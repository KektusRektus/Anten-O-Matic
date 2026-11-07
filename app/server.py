# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""Anten-O-Matic local GUI server: serves the UI on 127.0.0.1 and runs the evolution."""
import json, math, mimetypes, os, subprocess, sys, threading, time, traceback, webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs, unquote

APP = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(APP)                       # install folder; nothing is written outside it
OUT = os.path.join(ROOT, 'outputs')
DATA = os.path.join(ROOT, 'data')
DESIGNS = os.path.join(ROOT, 'designs')
os.environ.setdefault('MPLCONFIGDIR', os.path.join(DATA, 'matplotlib'))

import numpy as np
import core
import classic
import cad
import cad_evolve
import ems
MODULES = os.path.join(ROOT, 'modules')

STATE = {'running': False, 'progress': None, 'result': None, 'error': None, 'params': None, 'log': [], 'started': None}
LOCK = threading.Lock()
STOP = threading.Event()
LAST_SEEN = [time.time()]


def jsonable(o):
    if isinstance(o, dict): return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, (complex, np.complexfloating)): return [jsonable(o.real), jsonable(o.imag)]
    if isinstance(o, np.ndarray): return jsonable(o.tolist())
    return o


def log(msg):
    # live log for the window; the console copy is ASCII-safe (cp1252 consoles choke on Czech)
    with LOCK:
        STATE['log'].append([time.time(), msg])
        del STATE['log'][:-400]
        if STATE['progress'] is not None: STATE['progress'] = dict(STATE['progress'], beat=time.time())
    try:
        utf = 'utf' in (getattr(sys.stdout, 'encoding', '') or '').lower()
        print(time.strftime('%H:%M:%S'), msg if utf else msg.encode('ascii', 'replace').decode(), flush=True)
    except Exception: pass


LIVEPAT = {'fn': None, 'done': None, 'v': 0, 'pattern': None, 'ev': threading.Event()}


def progress(p):
    p = dict(p); fn = p.pop('_pat', None)
    with LOCK:
        STATE['progress'] = dict(p, beat=p.get('beat') or time.time(), livepat=LIVEPAT['v'])
        if fn is not None and fn is not LIVEPAT['done']: LIVEPAT['fn'] = fn; LIVEPAT['ev'].set()


def _livepat_worker():
    """3D pattern of the current best shape, computed beside the evolution (newest best wins,
    at most one every 1.5 s) so the 3D view can show it live."""
    while True:
        LIVEPAT['ev'].wait(); LIVEPAT['ev'].clear()
        with LOCK: fn = LIVEPAT['fn']; LIVEPAT['fn'] = None
        if fn is None: continue
        t = time.time()
        try: pat = jsonable(fn())
        except Exception: pat = None
        with LOCK:
            LIVEPAT['done'] = fn
            if pat:
                LIVEPAT['pattern'] = pat; LIVEPAT['v'] += 1
                if STATE['progress'] is not None: STATE['progress'] = dict(STATE['progress'], livepat=LIVEPAT['v'])
        time.sleep(max(0, 1.5 - (time.time() - t)))


threading.Thread(target=_livepat_worker, daemon=True).start()


RUN = {'id': 0, 'ev': None}   # current run; a run abandoned by the Stop watchdog must not write results later


def load_pattern(out):
    try:
        with open(os.path.join(out, 'pattern.json')) as f: return json.load(f)
    except Exception:
        return None


def save_result(ev, best):
    g, f, info = best
    with LOCK: STATE['progress'] = dict(STATE['progress'] or {}, phase='save', msg='Saving outputs and render', beat=time.time())
    log('Saving outputs: render, report, bending instructions' + (', STEP model' if ev.d.get('step', True) else ''))
    out, text, summ, steps = core.save_outputs(g, ev.d, OUT)
    log('Saved to %s' % out)
    rel = os.path.relpath(out, OUT).replace(os.sep, '/')
    with LOCK:
        STATE['result'] = {'pattern': load_pattern(out), 'dir': out, 'rel': rel, 'report': text, 'summary': summ, 'fit': f, 'info': info,
                           'plan': core.plan_text(ev.d, g), 'wires': core.wires_json(g, ev.d),
                           'warnings': core.warnings(ev.d), 'params': ev.d}


def worker(params, rid):
    try:
        ev = core.Evolver(params, progress, STOP, log)
        RUN['ev'] = ev
        best = ev.run()
        if RUN['id'] != rid: return          # the Stop watchdog already took over this run
        if best is None or best[0] is None:
            raise RuntimeError('Stopped before the first design existed.')
        save_result(ev, best)
    except Exception as e:
        if RUN['id'] != rid: return
        traceback.print_exc()
        log('ERROR: %s' % e)
        with LOCK: STATE['error'] = str(e)
    finally:
        if RUN['id'] == rid:
            with LOCK: STATE['running'] = False


def stop_watchdog(rid, wait=15.0):
    """Stop pressed: the evolution normally ends within a second or two. If it is still busy after
    `wait` seconds (stuck NEC process, hung pool...), kill the worker processes, save the best shape
    found so far and hand the window back to the user."""
    t = time.time()
    while time.time() - t < wait:
        time.sleep(0.25)
        with LOCK: busy = STATE['running']
        if not busy or RUN['id'] != rid: return
    ev = RUN['ev']
    log('Stop watchdog: the evolution has not responded for %.0f s, killing the worker processes' % wait)
    RUN['id'] += 1                            # orphan the stuck run
    try:
        import multiprocessing
        if ev is not None and getattr(ev, 'pool', None) is not None: ev.pool.terminate()
        for c in multiprocessing.active_children(): c.terminate()
    except Exception as e:
        log('Watchdog: killing the processes failed: %s' % e)
    try:
        best = ev.best if ev is not None else None
        if best and best[0] is not None:
            log('Saving the best shape found before the stop')
            if getattr(ev, 'kind', '') == 'cad': save_cad(best[0], evolved=True)
            else: save_result(ev, best)
        else:
            with LOCK: STATE['error'] = 'Stopped before the first design existed.'
    except Exception as e:
        traceback.print_exc()
        with LOCK: STATE['error'] = 'Stopped, but saving failed: %s' % e
    finally:
        with LOCK: STATE['running'] = False


def worker_classic(t, p):
    try:
        out_root = OUT
        if classic.TYPES[t][1] == 'calc':
            q = classic.prepare(t, p)
            log('Computing %s from formulas for %.3f MHz' % (classic.TYPES[t][0], q['f0']))
            progress({'phase': 'save', 'msg': 'Computing and drawing', 'gen': 1, 'gens': 1})
            out, text = (classic.patch_save if t == 'patch' else classic.dish_save)(q, out_root)
            res = {'kind': 'classic', 'type': t, 'calc': True, 'report': text, 'build': text, 'summary': [], 'wires': [], 'params': q}
        else:
            tu = classic.Tuner(t, p, progress, STOP, log)
            q = tu.run()
            with LOCK: STATE['progress'] = dict(STATE['progress'] or {}, phase='save', msg='Saving outputs and render', beat=time.time())
            log('Computing the band sweep, radiation patterns, render' + (' and STEP model' if q.get('step', True) else '') + '...')
            out, text, summ, btxt, wj = classic.save(t, q, out_root)
            res = {'kind': 'classic', 'type': t, 'calc': False, 'report': text, 'build': btxt, 'summary': summ, 'wires': wj,
                   'params': q, 'ground': t == 'helix'}
        log('Saved to %s' % out)
        res['pattern'] = load_pattern(out)
        res.update(dir=out, rel=os.path.relpath(out, OUT).replace(os.sep, '/'), warnings=[], label=classic.TYPES[t][0])
        with LOCK: STATE['result'] = res
    except Exception as e:
        traceback.print_exc()
        log('ERROR: %s' % e)
        with LOCK: STATE['error'] = str(e)
    finally:
        with LOCK: STATE['running'] = False


def ems_exe():
    exe = ems.find_exe(ROOT)
    if not exe:
        raise RuntimeError('openEMS not found. Download the official openEMS Windows zip and unpack it so that '
                           'openEMS.exe is in %s' % os.path.join(MODULES, 'openEMS'))
    return exe


def save_cad(doc, evolved=False):
    with LOCK: STATE['progress'] = dict(STATE['progress'] or {}, phase='save', msg='Saving outputs and render', beat=time.time())
    if doc.get('engine') == 'openems':
        log('Simulating in openEMS (FDTD): this takes minutes, the log shows the progress')
        st = lambda m, pct: progress(dict(STATE['progress'] or {}, phase='ems', msg=m, gen=round(pct, 1), gens=100, beat=time.time()))
        out, text, summ, btxt, wj, warn = ems.save(doc, OUT, ems_exe(), log, STOP, status=st)
    else:
        out, text, summ, btxt, wj, warn = cad.save(doc, OUT, log)
    log('Saved to %s' % out)
    keep = {k: v for k, v in doc.items() if not k.startswith('_')}
    res = {'kind': 'cad', 'type': 'cad', 'calc': False, 'report': text, 'build': btxt, 'summary': summ, 'wires': wj,
           'params': dict(keep, wire_d=0.001), 'ground': doc['ground'] == 'perfect', 'warnings': warn, 'doc': keep if evolved else None,
           'pattern': load_pattern(out), 'dir': out, 'rel': os.path.relpath(out, OUT).replace(os.sep, '/'), 'label': doc['name']}
    with LOCK: STATE['result'] = res


def worker_cad_evolve(doc, goal, opts, rid):
    try:
        ev = cad_evolve.CadEvolver(doc, goal, opts.get('pop', 40), opts.get('gens', 30), opts.get('jobs', 0), opts.get('seed', 0), progress, STOP, log)
        RUN['ev'] = ev
        best = ev.run()
        if RUN['id'] != rid: return
        if not best: raise RuntimeError('Stopped before the first design existed.')
        save_cad(best[0], evolved=True)
    except Exception as e:
        if RUN['id'] != rid: return
        traceback.print_exc()
        log('ERROR: %s' % e)
        with LOCK: STATE['error'] = str(e)
    finally:
        if RUN['id'] == rid:
            with LOCK: STATE['running'] = False


def worker_cad(doc):
    try:
        eng = 'openEMS (FDTD)' if doc.get('engine') == 'openems' else 'NEC-2'
        progress({'phase': 'save', 'msg': 'Simulating the CAD design in %s' % eng, 'gen': 1, 'gens': 1})
        log('Simulating %s in %s: %d elements, %.3f MHz' % (doc['name'], eng, len(doc['elements']), doc['f0']))
        save_cad(doc)
    except Exception as e:
        traceback.print_exc()
        log('ERROR: %s' % e)
        with LOCK: STATE['error'] = str(e)
    finally:
        with LOCK: STATE['running'] = False


def safe_name(n):
    n = ''.join(c if c.isalnum() or c in '-_ ' else '_' for c in str(n or '')).strip()[:60]
    if not n: raise ValueError('the design needs a name')
    return n


def defaults(q):
    d = {k: v[0] for k, v in q.items()}
    for k in ('rmax', 'hmax', 'maxlen'): d[k] = 0
    d['wire_d'] = float(d.get('wire_d', 2)) / 1000
    n = core.normalize(d)
    return {'rmax': n['rmax'] * 1000, 'hmax': n['hmax'] * 1000, 'maxlen': n['maxlen'] * 1000,
            'lam': core.lam(n) * 1000, 'warnings': core.warnings(n), 'jobs': n['jobs']}


class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass

    def send(self, code, body, ctype='application/json; charset=utf-8'):
        if not isinstance(body, (bytes, bytearray)):
            body = json.dumps(jsonable(body), ensure_ascii=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def file(self, base, rel):
        path = os.path.realpath(os.path.join(base, unquote(rel)))
        if not path.startswith(os.path.realpath(base) + os.sep) or not os.path.isfile(path):
            return self.send(404, {'error': 'not found'})
        ctype = mimetypes.guess_type(path)[0] or 'application/octet-stream'
        if path.endswith('.js'): ctype = 'text/javascript; charset=utf-8'
        with open(path, 'rb') as f: self.send(200, f.read(), ctype)

    def do_GET(self):
        LAST_SEEN[0] = time.time()
        u = urlparse(self.path)
        if u.path == '/': return self.file(os.path.join(APP, 'static'), 'index.html')
        if u.path.startswith('/static/'): return self.file(os.path.join(APP, 'static'), u.path[8:])
        if u.path.startswith('/out/'): return self.file(OUT, u.path[5:])
        if u.path == '/api/status':
            with LOCK: return self.send(200, dict(STATE, now=time.time()))
        if u.path == '/api/defaults':
            try: return self.send(200, defaults(parse_qs(u.query)))
            except Exception as e: return self.send(400, {'error': str(e)})
        if u.path == '/api/livepat':
            with LOCK: return self.send(200, {'v': LIVEPAT['v'], 'pattern': LIVEPAT['pattern']})
        if u.path == '/api/ping': return self.send(200, {'ok': True})
        if u.path == '/api/materials': return self.send(200, cad.load_materials(DATA))
        if u.path == '/api/ems': return self.send(200, {'exe': ems.find_exe(ROOT), 'folder': os.path.join(MODULES, 'openEMS')})
        if u.path == '/api/designs':
            os.makedirs(DESIGNS, exist_ok=True)
            return self.send(200, sorted(f[:-5] for f in os.listdir(DESIGNS) if f.endswith('.json')))
        if u.path == '/api/design':
            try:
                n = safe_name(parse_qs(u.query).get('name', [''])[0])
                with open(os.path.join(DESIGNS, n + '.json'), encoding='utf-8') as f: return self.send(200, json.load(f))
            except Exception as e: return self.send(404, {'error': 'cannot load the design: %s' % e})
        if u.path == '/api/classic_types':
            q = parse_qs(u.query)
            try: return self.send(200, classic.catalogue(float(q.get('f0', ['145'])[0]), float(q.get('wire_d', ['2'])[0])))
            except Exception as e: return self.send(400, {'error': str(e)})
        self.send(404, {'error': 'not found'})

    def do_POST(self):
        LAST_SEEN[0] = time.time()
        u = urlparse(self.path)
        n = int(self.headers.get('Content-Length') or 0)
        body = json.loads(self.rfile.read(n) or b'{}')
        if u.path == '/api/ems_mesh':      # mesh size estimate for the Simulation panel
            try: return self.send(200, ems.mesh_info(cad.clean_doc(dict(body.get('doc') or {}, engine='openems'), cad.load_materials(DATA))))
            except Exception as e: return self.send(200, {'error': str(e)})
        if u.path == '/api/materials':
            try: return self.send(200, cad.save_materials(DATA, body.get('materials') or []))
            except Exception as e: return self.send(400, {'error': str(e)})
        if u.path == '/api/design':
            try:
                n = safe_name(body.get('name'))
                os.makedirs(DESIGNS, exist_ok=True)
                with open(os.path.join(DESIGNS, n + '.json'), 'w', encoding='utf-8') as f: json.dump(dict(body.get('doc') or {}, name=n), f, indent=1, ensure_ascii=False)
                return self.send(200, {'ok': True, 'name': n})
            except Exception as e: return self.send(400, {'error': str(e)})
        if u.path == '/api/pcb_svg':            # planar layout -> outputs/svg/<name>_<time>.svg (never outside the outputs folder)
            try:
                e = dict(body.get('element') or {}); e['type'] = 'pcb'; e.setdefault('id', 'svg'); e.setdefault('mat', 'Copper')
                doc = cad.clean_doc({'name': 'svg', 'f0': 1000, 'bw': 0, 'elements': [e]}, cad.load_materials(DATA))
                el = doc['elements'][0]
                d = os.path.join(OUT, 'svg'); os.makedirs(d, exist_ok=True)
                fn = os.path.join(d, '%s_%s.svg' % (safe_name(el['name']), time.strftime('%Y%m%d-%H%M%S')))
                with open(fn, 'w', encoding='utf-8') as f: f.write(cad.PCB.svg(el))
                return self.send(200, {'path': fn})
            except Exception as e: return self.send(400, {'error': str(e)})
        if u.path == '/api/classic_to_cad':
            try:
                t = body.get('type')
                if t not in classic.TYPES or classic.TYPES[t][1] == 'calc': raise ValueError('only wire antennas can be opened in the CAD')
                p = dict(body.get('params') or {})
                if float(p.get('wire_d') or 0) > 0.5: p['wire_d'] = float(p['wire_d']) / 1000   # mm -> m
                p = classic.prepare(t, {k: v for k, v in p.items() if k != '_lam'})
                doc, notes = classic.to_cad(t, p)
                return self.send(200, {'doc': doc, 'notes': notes})
            except Exception as e:
                traceback.print_exc()
                return self.send(400, {'error': str(e)})
        if u.path == '/api/run' and body.get('kind') == 'cad':
            with LOCK:
                if STATE['running']: return self.send(409, {'error': 'Already running.'})
                try:
                    doc = cad.clean_doc(body.get('doc'), cad.load_materials(DATA))
                    if doc['engine'] == 'openems':
                        ems_exe(); ems.geometry(doc)
                    else:
                        cad.nec_wires(doc)
                    evo = body.get('evolve')
                    if evo:
                        goal = {'kind': evo.get('goal', 'match'), 'az': float(evo.get('az', 0)), 'el': float(evo.get('el', 0)), 'pol': evo.get('pol', 'total')}
                        if goal['kind'] not in ('match', 'dir', 'target'): raise ValueError('unknown goal')
                        if goal['pol'] not in ('total', 'rhcp', 'lhcp'): raise ValueError('unknown polarisation')
                        if goal['kind'] == 'target': goal['target'] = core.PAT.clean_target(evo.get('target'))
                        if doc['engine'] == 'openems':
                            goal['exe'] = ems_exe(); goal['simroot'] = os.path.join(OUT, '_evolve_openems')
                        if cad_evolve.Space(doc).n == 0: raise ValueError('nothing may change: give at least one element evolution freedom (shape, free-form or pixels), or a dimension a range')
                        opts = {k: int(evo.get(k) or 0) for k in ('pop', 'gens', 'jobs', 'seed')}
                        opts['pop'] = opts['pop'] or 40; opts['gens'] = opts['gens'] or 30
                except Exception as e:
                    return self.send(400, {'error': 'Invalid design: %s' % e})
                STOP.clear()
                STATE.update(running=True, progress=None, result=None, error=None, params=body, log=[], started=time.time())
            RUN['id'] += 1; RUN['ev'] = None
            if evo: threading.Thread(target=worker_cad_evolve, args=(doc, goal, opts, RUN['id']), daemon=True).start()
            else: threading.Thread(target=worker_cad, args=(doc,), daemon=True).start()
            return self.send(200, {'ok': True})
        if u.path == '/api/run':
            if body.get('kind') == 'classic':
                with LOCK:
                    if STATE['running']: return self.send(409, {'error': 'Already running.'})
                    try:
                        t = body['type']
                        if t not in classic.TYPES: raise ValueError('unknown type %s' % t)
                        p = dict(body.get('params') or {})
                        for k in ('name', 'f0', 'bw', 'z0', 'swr', 'material', 'step', 'iters', 'jobs', 'pol', 'seed', 'bend_r'):
                            if k in body: p[k] = body[k]
                        p['wire_d'] = float(body.get('wire_d') or 2) / 1000
                        if p.get('f0') in (None, '') or float(p['f0']) <= 0: raise ValueError('frequency must be positive')
                        classic.prepare(t, p)
                    except Exception as e:
                        return self.send(400, {'error': 'Invalid parameter: %s' % e})
                    STOP.clear()
                    STATE.update(running=True, progress=None, result=None, error=None, params=body, log=[], started=time.time())
                RUN['id'] += 1; RUN['ev'] = None
                threading.Thread(target=worker_classic, args=(t, p), daemon=True).start()
                return self.send(200, {'ok': True})
            with LOCK:
                if STATE['running']: return self.send(409, {'error': 'Already running.'})
                try:
                    p = dict(body)
                    for k in ('rmax', 'hmax', 'maxlen', 'wire_d', 'bend_r'):
                        p[k] = float(p.get(k) or 0) / 1000
                    core.normalize(p)
                except Exception as e:
                    return self.send(400, {'error': 'Invalid parameter: %s' % e})
                STOP.clear()
                STATE.update(running=True, progress=None, result=None, error=None, params=body, log=[], started=time.time())
            RUN['id'] += 1; RUN['ev'] = None
            threading.Thread(target=worker, args=(p, RUN['id']), daemon=True).start()
            return self.send(200, {'ok': True})
        if u.path == '/api/stop':
            first = not STOP.is_set()
            STOP.set()
            with LOCK: busy = STATE['running']
            if busy and first:
                log('Stop received, finishing the work in progress')
                threading.Thread(target=stop_watchdog, args=(RUN['id'],), daemon=True).start()
            return self.send(200, {'ok': True, 'running': busy})
        if u.path == '/api/open':
            with LOCK: r = STATE['result']
            target = r['dir'] if r else OUT
            os.makedirs(target, exist_ok=True)
            try:
                if sys.platform == 'win32': os.startfile(target)
                elif sys.platform == 'darwin': subprocess.Popen(['open', target])
                else: subprocess.Popen(['xdg-open', target])
            except Exception as e:
                return self.send(500, {'error': str(e)})
            return self.send(200, {'ok': True})
        if u.path == '/api/quit':
            self.send(200, {'ok': True}); threading.Timer(0.3, lambda: os._exit(0)).start(); return
        self.send(404, {'error': 'not found'})


def find_browser():
    if sys.platform != 'win32': return None
    cands = []
    for env in ('PROGRAMFILES(X86)', 'PROGRAMFILES', 'LOCALAPPDATA'):
        base = os.environ.get(env)
        if not base: continue
        cands += [os.path.join(base, 'Microsoft', 'Edge', 'Application', 'msedge.exe'),
                  os.path.join(base, 'Google', 'Chrome', 'Application', 'chrome.exe')]
    return next((c for c in cands if os.path.isfile(c)), None)


def open_window(url):
    b = find_browser()
    if b:
        prof = os.path.join(DATA, 'window')  # browser profile kept inside the install folder
        os.makedirs(prof, exist_ok=True)
        subprocess.Popen([b, '--app=' + url, '--user-data-dir=' + prof, '--window-size=1440,920',
                          '--no-first-run', '--no-default-browser-check'])
    else:
        webbrowser.open(url)


def watchdog():
    # quit when the window has been closed (no polling for a while) and nothing is running
    while True:
        time.sleep(5)
        with LOCK: busy = STATE['running']
        if not busy and time.time() - LAST_SEEN[0] > 60:
            os._exit(0)


def main(open_ui=True, port=0):
    for d_ in (OUT, DATA, DESIGNS, MODULES): os.makedirs(d_, exist_ok=True)
    if sys.stdout is None or sys.stderr is None:  # pythonw: keep a log inside the install folder
        log = open(os.path.join(DATA, 'log.txt'), 'a', encoding='utf-8', buffering=1)
        sys.stdout = sys.stderr = log
    srv = ThreadingHTTPServer(('127.0.0.1', port), H)
    url = 'http://127.0.0.1:%d/' % srv.server_address[1]
    print('Anten-O-Matic running at', url, flush=True)
    if open_ui:
        open_window(url)
        threading.Thread(target=watchdog, daemon=True).start()
    srv.serve_forever()


if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()
    main(open_ui='--no-window' not in sys.argv,
         port=int(sys.argv[sys.argv.index('--port') + 1]) if '--port' in sys.argv else 0)
