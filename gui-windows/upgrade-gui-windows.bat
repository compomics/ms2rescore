@echo off
call Miniconda3/condabin/conda.bat activate ms2rescore
echo Updating MS²Rescore...
call pip install ms2rescore[gui] --upgrade
echo.
echo Everything updated!
echo.
pause
