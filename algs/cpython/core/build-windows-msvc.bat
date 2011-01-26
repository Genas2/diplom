@del /Q alglib.dll > NUL 2> NUL
@copy .\src\*.c tmp > NUL 2> NUL
@copy .\src\*.h  tmp > NUL 2> NUL
@cd tmp
@cl /LD /Fealglib.dll /DX_FOR_WINDOWS /DAE_USE_ALLOC_COUNTER /I. %~1 *.c 
@if not ERRORLEVEL 1 goto lbl_build_ok
@cd ..
@exit /B 1
:lbl_build_ok
@copy alglib.dll ..\alglib.dll > NUL 2> NUL
@del /Q *
@cd ..
