@ECHO OFF
TITLE Finding Liutex Core Lines
gfortran -c .\liutex_mod.f03
gfortran -c .\coreline_v1.f03
gfortran .\coreline_v1.o .\liutex_mod.o -o find_core_lines.exe
START /WAIT /B "Finding Liutex Core Lines" ".\find_core_lines.exe"
PAUSE