@echo off
chcp 65001 >nul
set MATLAB="C:\Program Files\MATLAB\R2025a\bin\matlab.exe"
set WORKDIR=c:\Users\z5242\Desktop\work\Nonlinear_Systems-main\Nonlinear_Systems-main\1_Lyapunov_method_YALMIP\新建文件夹
%MATLAB% -batch "cd('%WORKDIR%'); work24_6dof" > "%WORKDIR%\matlab_log.txt" 2>&1
echo EXIT_CODE=%ERRORLEVEL% >> "%WORKDIR%\matlab_log.txt"
