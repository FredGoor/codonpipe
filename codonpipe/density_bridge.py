"""codonpipe.density_bridge

Utility to execute Plotting_Pipeline.py programmatically.
"""

# codonpipe/density_bridge.py
import os
import sys
import uuid
import importlib.util
import traceback


def run_density_plot_script(excel_path, pipeline_cfg):
    """
    Executes the plotting pipeline script (e.g. Plotting_Pipeline.py) as a standalone module.

    IMPORTANT DESIGN CHOICE:
    - This bridge ONLY wires the workbook + required sheet names.
    - It does NOT override any plotting aesthetics or COLOR_MODE/enrichment settings.
      Those must live in the plotting script itself (Plotting_Pipeline.py).
    """
    script_path = pipeline_cfg.get("plotting_pipeline_script_path", "")
    if not script_path:
        print("[WARN] No density_plot_script_path specified; skipping plotting pipeline.")
        return

    script_path = os.path.abspath(os.path.expanduser(script_path))
    if not os.path.exists(script_path):
        print(f"[WARN] Plotting pipeline script not found at:\n  {script_path}\nSkipping.")
        return

    script_name = os.path.basename(script_path)
    print(f"[INFO] Running plotting pipeline: {script_name}")

    # Always load as a fresh module to avoid interactive environments (e.g. IPython) caching effects
    module_name = f"_codonpipe_plot_{uuid.uuid4().hex}"
    try:
        spec = importlib.util.spec_from_file_location(module_name, script_path)
        if spec is None or spec.loader is None:
            raise RuntimeError("Could not create import spec for plotting script.")

        mod = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = mod
        spec.loader.exec_module(mod)

        # -------- Required wiring (only) --------
        if hasattr(mod, "EXCEL_PATH"):
            mod.EXCEL_PATH = excel_path
        if hasattr(mod, "USE_FILE_DIALOG"):
            mod.USE_FILE_DIALOG = False
        if hasattr(mod, "UMAP_SHEET") and "sheet_coordinates" in pipeline_cfg:
            mod.UMAP_SHEET = pipeline_cfg["sheet_coordinates"]
        if hasattr(mod, "DATA_SHEET") and "sheet_binary" in pipeline_cfg:
            mod.DATA_SHEET = pipeline_cfg["sheet_binary"]

        # -------- Debug prints (truth source) --------
        print(f"[INFO]   EXCEL_PATH : {getattr(mod, 'EXCEL_PATH', '(missing)')}")
        print(f"[INFO]   UMAP_SHEET : {getattr(mod, 'UMAP_SHEET', '(missing)')}")
        print(f"[INFO]   DATA_SHEET : {getattr(mod, 'DATA_SHEET', '(missing)')}")
        if hasattr(mod, "COLOR_MODE"):
            print(f"[INFO]   {script_name}.COLOR_MODE = {getattr(mod, 'COLOR_MODE')}")
        else:
            print(f"[INFO]   {script_name} has no global COLOR_MODE (may be inside cfg only).")
        print("[INFO]   NOTE: COLOR_MODE/enrichment settings are NOT overridden by the bridge.")

        # -------- Execute --------
        if hasattr(mod, "main") and callable(mod.main):
            mod.main()
        else:
            raise RuntimeError(f"{script_name} has no callable main(); nothing executed.")

    except Exception as e:
        print(f"[WARN] Failed to run plotting pipeline script: {e}")
        print(traceback.format_exc())
