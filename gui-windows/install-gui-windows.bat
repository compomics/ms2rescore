@echo off

set downloadurl=https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
echo Downloading Miniconda...
powershell -Command "Invoke-WebRequest %downloadurl% -OutFile Miniconda3-latest-Windows-x86_64.exe"

echo Installing Miniconda...
pushd %~dp0
set script_dir=%CD%
start /wait "" Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /S /D=%script_dir%\Miniconda3
echo Installing new environment in Miniconda...
call Miniconda3/Scripts/activate.bat
call conda create --yes -n ms2rescore python=3.8
call Miniconda3/Scripts/activate.bat ms2rescore
call conda config --add channels defaults
call conda config --add channels bioconda
call conda config --add channels conda-forge
echo Installing MS²Rescore...
call pip install ms2rescore[gui]
echo.
echo Everything installed!
echo.
echo Do not forget to install Percolator and percolator-convertors seperately
echo https://github.com/percolator/percolator/releases/latest
echo.
pause
