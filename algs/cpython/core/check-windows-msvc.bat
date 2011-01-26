@copy .\src\*.c tmp > NUL 2> NUL
@copy .\src\*.h  tmp > NUL 2> NUL
@copy .\tests\*.c tmp > NUL 2> NUL
@copy .\tests\*.h  tmp > NUL 2> NUL
@cd tmp
@cl /Fecheck-alglib.exe /DX_FOR_WINDOWS /DAE_USE_ALLOC_COUNTER /I. %~1 *.c 
@if not ERRORLEVEL 1 goto lbl_build_ok
@cd ..
@exit /B 1
:lbl_build_ok
check-alglib.exe
@if not ERRORLEVEL 1 goto lbl_check_ok
@echo Errors detected
@cd ..
@exit /B 1
:lbl_check_ok
@del /Q *
@cd ..
