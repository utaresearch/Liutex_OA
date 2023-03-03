@ECHO OFF
TITLE Finding Liutex Core Line
gfortran -c .\liutex_mod.f03
gfortran -c .\coreline_v1.f03
gfortran .\coreline_v1.o .\liutex_mod.o -o find_core_lines.exe
START "" ".\find_core_lines.exe"
EXIT