# Anten-O-Matic. Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors
# SPDX-License-Identifier: GPL-3.0-or-later (see LICENSE.txt)
"""ctypes binding to the necpp C library (necpp.dll on Windows, libnecpp.so elsewhere)."""
import ctypes, os, sys

_here = os.path.dirname(os.path.abspath(__file__))
_name = 'necpp.dll' if sys.platform == 'win32' else 'libnecpp.so'
_lib = ctypes.CDLL(os.path.join(_here, _name))

P, I, D = ctypes.c_void_p, ctypes.c_int, ctypes.c_double
_sig = {
    'nec_create': (P, []),
    'nec_delete': (ctypes.c_long, [P]),
    'nec_wire': (ctypes.c_long, [P, I, I, D, D, D, D, D, D, D, D, D]),
    'nec_geometry_complete': (ctypes.c_long, [P, I]),
    'nec_gn_card': (ctypes.c_long, [P, I, I, D, D, D, D, D, D]),
    'nec_fr_card': (ctypes.c_long, [P, I, I, D, D]),
    'nec_ld_card': (ctypes.c_long, [P, I, I, I, I, D, D, D]),
    'nec_excitation_voltage': (ctypes.c_long, [P, I, I, D, D]),
    'nec_tl_card': (ctypes.c_long, [P, I, I, I, I, D, D, D, D, D, D]),
    'nec_rp_card': (ctypes.c_long, [P, I, I, I, I, I, I, I, D, D, D, D, D, D]),
    'nec_xq_card': (ctypes.c_long, [P, I]),
    'nec_ek_card': (ctypes.c_long, [P, I]),
    'nec_gain': (D, [P, I, I, I]),
    'nec_gain_point_rhcp': (D, [P, I, I, I]),
    'nec_gain_point_lhcp': (D, [P, I, I, I]),
    'nec_error_message': (ctypes.c_char_p, []),
    'nec_impedance_real_indexed': (D, [P, I, I]),
    'nec_impedance_imag_indexed': (D, [P, I, I]),
}
for pol in ('', 'rhcp_', 'lhcp_'):
    for st in ('max', 'min', 'mean'):
        _sig['nec_gain_%s%s' % (pol, st)] = (D, [P, I])

for fn, (res, args) in _sig.items():
    f = getattr(_lib, fn)
    f.restype, f.argtypes = res, args
    globals()[fn] = f


def gain(n, pol, stat, fi):
    """pol: 'rhcp' | 'lhcp' | 'total'; stat: 'max' | 'min' | 'mean'."""
    name = 'nec_gain_%s%s' % ('' if pol == 'total' else pol + '_', stat)
    return globals()[name](n, fi)
