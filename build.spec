import gooey


def get_version(path):
    """Get __version__ from Python file."""
    with open(path, "rt") as f:
        for line in f:
            if line.startswith("__version__ = "):
                return line.strip().split(" = ")[1].strip("\"'")

# Gooey
gooey_root = os.path.dirname(gooey.__file__)
gooey_languages = Tree(os.path.join(gooey_root, 'languages'), prefix = 'gooey/languages')
gooey_images = Tree(os.path.join(gooey_root, 'images'), prefix = 'gooey/images')
image_overrides = Tree('ms2rescore/package_data/img', prefix='ms2rescore/package_data/img')

# Package data and version
package_data = Tree('ms2rescore/package_data', prefix='ms2rescore/package_data')
version = get_version("ms2rescore/_version.py")

a = Analysis(['ms2rescore/gui.py'],
             pathex=['c:\\Python27\\Scripts'],
             hiddenimports=["ms2pip.ms2pipC", "xgboost"],
             hookspath="pyinstaller/hooks",
             runtime_hooks=None,
             )
pyz = PYZ(a.pure)

options = [('u', None, 'OPTION'), ('u', None, 'OPTION'), ('u', None, 'OPTION')]

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          options,
          gooey_languages,
          gooey_images,
          image_overrides,
          package_data,
          name='ms2rescore_' + version,
          debug=False,
          strip=None,
          upx=True,
          console=False,
          windowed=True,
          icon='ms2rescore/package_data/img/program_icon.ico')
