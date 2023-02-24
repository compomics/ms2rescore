import os
import sys

import importlib_metadata
from PyInstaller.utils.hooks import collect_all
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, Tree
from ms2rescore import __version__

# Package info
exe_name = "ms2rescore"
script_name = "ms2rescore/gui.py"
icon = './img/ms2rescore.ico'
location = os.getcwd()
project = "ms2rescore"
bundle_name = "ms2rescore"
bundle_identifier = f"{bundle_name}.{__version__}"


# Collect hidden imports and datas for all requirements
requirements = {req.split()[0] for req in importlib_metadata.requires(project)}
requirements.update([project, "distributed", "customtkinter", "ms2pip"])
hidden_imports = set()
datas = []
binaries = []
checked = set()
while requirements:
    requirement = requirements.pop()
    checked.add(requirement)
    try:
        module_version = importlib_metadata.version(requirement)
    except (importlib_metadata.PackageNotFoundError, ModuleNotFoundError, ImportError):
        if requirement != "sklearn":
            continue
    try:
        datas_, binaries_, hidden_imports_ = collect_all(
            requirement, include_py_files=True
        )
    except ImportError:
        continue
    datas += datas_
    hidden_imports_ = set(hidden_imports_)
    if "" in hidden_imports_:
        hidden_imports_.remove("")
    if None in hidden_imports_:
        hidden_imports_.remove(None)
    requirements |= hidden_imports_ - checked
    hidden_imports |= hidden_imports_

hidden_imports = sorted([h for h in hidden_imports if "tests" not in h.split(".")])
hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]
datas = [d for d in datas if ("__pycache__" not in d[0]) and (d[1] not in [".", "build", "dist", "Output"])]
datas += [("ms2rescore\package_data", "package_data" )]

block_cipher = None
# Build package
a = Analysis(
    [script_name],
    pathex=[location],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False
)

pyz = PYZ(
    a.pure,
    a.zipped_data,
    cipher=block_cipher
)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name=exe_name,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    windowed=True,
    disable_windowed_traceback=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="./img/ms2rescore.ico"
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name=exe_name
)