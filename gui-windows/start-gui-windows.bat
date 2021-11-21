@echo off
call Miniconda3/condabin/conda.bat activate ms2rescore
start /b python -m ms2rescore.gui
