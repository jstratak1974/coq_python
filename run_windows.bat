@echo off
cd /d "%~dp0"
py coq_python.py
if errorlevel 1 pause
