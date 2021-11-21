@echo off
call Miniconda3/Scripts/activate.bat ms2rescore
start /B python -m ms2rescore.gui
