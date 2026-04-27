# -*- mode: python ; coding: utf-8 -*-
"""\
PyInstaller Spec File for OmicsVisStat
=========================================================

Build:
    pyinstaller statistics.spec

Output (one-folder build):
    dist/<app_name>/<app_name>.exe  # generated folder/executable name includes version

Runtime Requirements (Optional):
    - R (4.5.x or compatible) with limma package for covariate adjustment
      Install: In R console run:
        install.packages('BiocManager')
        BiocManager::install('limma')
    - rpy2 Python package (included if installed during build)
      Note: Linear Model (OLS) method works without R/rpy2
"""

from __future__ import annotations

import os
from PyInstaller.utils.hooks import collect_data_files, collect_submodules


block_cipher = None


def _safe_collect_data_files(package_name: str):
    try:
        return collect_data_files(package_name)
    except Exception:
        return []


def _safe_collect_submodules(package_name: str):
    try:
        return collect_submodules(package_name)
    except Exception:
        return []


# Get the current directory
spec_root = os.path.abspath(SPECPATH)

# version string for executable name (update manually or derive from release tag)
version = "1.0.0"
app_name = f"OmicsVisStat_v{version}"


# Data files
datas = []

# Core scientific/data stack (safe-guarded in case something isn't installed)
for pkg in [
    "numpy",
    "pandas",
    "scipy",
    "sklearn",
    "statsmodels",
    "matplotlib",
    "seaborn",
    "openpyxl",
    "xlsxwriter",
    # Optional
    "pyarrow",
    "plotly",
    "kaleido",
]:
    datas += _safe_collect_data_files(pkg)

# Bundle config assets (keep folder structure)
config_files = [
    ("main_script/factor_mapping_config.json", "main_script"),
    ("gui/tabs/statistics_config.json", "gui/tabs"),
    ("gui/tabs/visualization_config.json", "gui/tabs"),
    ("gui/tabs/ml_config.json", "gui/tabs"),
]

for src, dest in config_files:
    full_path = os.path.join(spec_root, src)
    if os.path.exists(full_path):
        datas.append((src, dest))


# Hidden imports
hiddenimports = [
    # Tkinter UI
    "tkinter",
    "tkinter.ttk",
    "tkinter.filedialog",
    "tkinter.messagebox",
    "tkinter.scrolledtext",
    "tkinter.simpledialog",

    # App modules
    "gui",
    "gui.main",
    "gui.shared",
    "gui.shared.base_tab",
    "gui.shared.column_assignment",
    "gui.shared.data_manager",
    "gui.shared.pairwise_column_mapper",
    "gui.shared.utils",
    "gui.tabs",
    "gui.tabs.statistics_tab",
    "gui.tabs.visualization_tab",
    "gui.tabs.machine_learning_tab",
    "gui.tabs.covariate_adjustment_section",
    "gui.tabs.help_tab",

    # main_script modules (directory is treated as a namespace package)
    "main_script.metabolite_statistics_analysis",
    "main_script.metabolites_visualization",
    "main_script.covariate_adjustment",
    "main_script.factor_mapping_manager",
    "main_script.one_way_anova",
    "main_script.nonparametric_two_way_anova",
    "main_script.two_way_anova_new_format",
    "main_script.two_way_anova_output_formatter",
    "main_script.help",
    
    # Encodings (CRITICAL: prevents startup crashes on different system locales)
    "encodings",
    "encodings.utf_8",
    "encodings.cp1252",
    "encodings.ascii",
    "encodings.latin_1",
    
    # NumPy core (often missed)
    "numpy.core._methods",
    "numpy.lib.format",
    "numpy.fft",
    "numpy.random",
    "numpy.linalg",
    
    # Pandas Excel support
    "pandas.io.excel._openpyxl",
    "pandas.io.excel._xlsxwriter",
    
    # PyArrow for feather files
    "pyarrow",
    "pyarrow.feather",
    
    # Platform specific
    "ctypes",
    "ctypes.wintypes",
]

# Pull in submodules for key dependencies (safe)
hiddenimports += _safe_collect_submodules("scipy")
hiddenimports += _safe_collect_submodules("sklearn")
hiddenimports += _safe_collect_submodules("statsmodels")
hiddenimports += _safe_collect_submodules("matplotlib")
hiddenimports += _safe_collect_submodules("numpy")
hiddenimports += _safe_collect_submodules("pandas")

# Collect binaries from packages with compiled extensions (CRITICAL for other machines)
from PyInstaller.utils.hooks import collect_dynamic_libs
binaries = []
binaries += collect_dynamic_libs('sklearn')
binaries += collect_dynamic_libs('scipy')
binaries += collect_dynamic_libs('numpy')


a = Analysis(
    ["run_gui.py"],
    pathex=[spec_root],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "IPython",
        "jupyter",
        "notebook",
        "PyQt5",
        "PyQt6",
        "PySide2",
        "PySide6",
        "wx",
        "pytest",
        "sphinx",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)


pyz = PYZ(
    a.pure,
    a.zipped_data,
    cipher=block_cipher,
)


exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="OmicsVisStat",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,  # DISABLED: UPX can corrupt scipy/numpy DLLs causing crashes on other machines
    console=True,  # ENABLED: Shows error messages when app crashes (set to False for release)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,
    version=None,
)


coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,  # DISABLED: UPX can corrupt scipy/numpy DLLs causing crashes on other machines
    upx_exclude=[],
    name="OmicsVisStat",
)

