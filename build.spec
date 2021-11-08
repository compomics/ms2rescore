import gooey

gooey_root = os.path.dirname(gooey.__file__)
gooey_languages = Tree(os.path.join(gooey_root, 'languages'), prefix = 'gooey/languages')
gooey_images = Tree(os.path.join(gooey_root, 'images'), prefix = 'gooey/images')
image_overrides = Tree('ms2rescore/package_data/img', prefix='ms2rescore/package_data/img')
package_data = Tree('ms2rescore/package_data', prefix='ms2rescore/package_data')
a = Analysis(['ms2rescore/gui.py'],
             pathex=['c:\\Python27\\Scripts'],
             hiddenimports=[],
             hookspath=None,
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
          gooey_languages, # Add them in to collected files
          gooey_images, # Same here.
          image_overrides,
          package_data,
          name='MS²Rescore',
          debug=False,
          strip=None,
          upx=True,
          console=False,
          windowed=True,
          icon='ms2rescore/package_data/img/program_icon.ico')
