Unicode true
!include "MUI2.nsh"
!define APP "Anten-O-Matic"
!define VER "3.7.2"
Name "${APP} ${VER}"
OutFile "Anten-O-Matic-Setup-${VER}.exe"
InstallDir "$PROFILE\Desktop\AITools\Anten-O-Matic"
RequestExecutionLevel user
SetCompressor /SOLID lzma
SetCompressorDictSize 64
!define MUI_ICON "anten-o-matic.ico"
!define MUI_UNICON "anten-o-matic.ico"
!define MUI_FINISHPAGE_RUN
!define MUI_FINISHPAGE_RUN_TEXT "Launch Anten-O-Matic"
!define MUI_FINISHPAGE_RUN_FUNCTION LaunchApp
!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\README.txt"
!define MUI_FINISHPAGE_SHOWREADME_NOTCHECKED
!define MUI_WELCOMEPAGE_TEXT "Evolved wire antennas: NEC-2 + a genetic algorithm.$\r$\n$\r$\nEverything (program and outputs) is installed into one folder. Nothing is written elsewhere and nothing is sent to the internet."
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_LANGUAGE "English"

Function LaunchApp
  ExecShell "" "$INSTDIR\Anten-O-Matic.lnk"
FunctionEnd

Section "Anten-O-Matic (required)" SecMain
  SectionIn RO
  SetOutPath "$INSTDIR"
  ; upgrade from 2.x: drop the old app and the Czech files, keep the user's outputs
  RMDir /r "$INSTDIR\app"
  Delete "$INSTDIR\CTIMNE.txt"
  Delete "$INSTDIR\Spustit.cmd"
  Delete "$INSTDIR\Odinstalovat.exe"
  IfFileExists "$INSTDIR\outputs\*.*" +3 0
  IfFileExists "$INSTDIR\vystupy\*.*" 0 +2
  Rename "$INSTDIR\vystupy" "$INSTDIR\outputs"
  File /r "python"
  File /r "app"
  File /r "licence"
  File "README.txt"
  File "LICENSE.txt"
  File "Start.cmd"
  File "anten-o-matic.ico"
  CreateDirectory "$INSTDIR\outputs"
  CreateDirectory "$INSTDIR\data"
  CreateDirectory "$INSTDIR\designs"
  CreateDirectory "$INSTDIR\modules"
  SetOutPath "$INSTDIR\modules"
  File "modules\README.txt"
  SetOutPath "$INSTDIR"
  CreateShortCut "$INSTDIR\Anten-O-Matic.lnk" "$INSTDIR\python\pythonw.exe" '"$INSTDIR\app\server.py"' "$INSTDIR\anten-o-matic.ico" 0 SW_SHOWNORMAL "" "Anten-O-Matic"
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section /o "Desktop shortcut (outside the install folder)" SecDesk
  CreateShortCut "$DESKTOP\Anten-O-Matic.lnk" "$INSTDIR\python\pythonw.exe" '"$INSTDIR\app\server.py"' "$INSTDIR\anten-o-matic.ico" 0
SectionEnd

Section "Uninstall"
  Delete "$DESKTOP\Anten-O-Matic.lnk"
  RMDir /r "$INSTDIR\python"
  RMDir /r "$INSTDIR\app"
  RMDir /r "$INSTDIR\licence"
  ; keep the user's material library (data\materials.json), designs, outputs and modules
  RMDir /r "$INSTDIR\data\window"
  RMDir /r "$INSTDIR\data\matplotlib"
  Delete "$INSTDIR\data\log.txt"
  RMDir "$INSTDIR\data"
  Delete "$INSTDIR\modules\README.txt"
  Delete "$INSTDIR\README.txt"
  Delete "$INSTDIR\LICENSE.txt"
  Delete "$INSTDIR\Start.cmd"
  Delete "$INSTDIR\anten-o-matic.ico"
  Delete "$INSTDIR\Anten-O-Matic.lnk"
  Delete "$INSTDIR\Uninstall.exe"
  RMDir "$INSTDIR"
SectionEnd
