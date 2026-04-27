#!/usr/bin/env python3
"""
Entry point for the OmicsVisStat GUI application.

This script:
1. Properly configures Python paths
2. Suppresses background thread errors
3. Launches the OmicsVisStat GUI
4. Handles cleanup on exit
"""
import sys
import os
import threading
import warnings

# Suppress thread-related warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Determine the correct base directory for PyInstaller or normal execution
if getattr(sys, 'frozen', False):
    # Running as compiled PyInstaller executable
    # Use the directory where the exe is located
    application_path = os.path.dirname(sys.executable)
else:
    # Running as normal Python script
    application_path = os.path.dirname(os.path.abspath(__file__))

# Change working directory to the application path
# This ensures relative paths work correctly when double-clicking the exe
os.chdir(application_path)

# Add the script directory to Python path
sys.path.insert(0, application_path)

# Note: Database files are not required for OmicsVisStat functionality
# The application can work with user-provided Excel/CSV files directly
print(f"Application path: {application_path}")
print(f"Current working directory: {os.getcwd()}")
print("Starting OmicsVisStat: Statistical and Visualization Platform for Omics Data Analysis...")

# Monkey-patch to suppress background thread errors
original_excepthook = sys.excepthook
def suppress_thread_errors(exc_type, exc_value, traceback):
    """Suppress 'main thread is not in main loop' errors from background threads"""
    # Suppress specific background thread error
    if (exc_type == RuntimeError and 
        "main thread is not in main loop" in str(exc_value) and
        threading.current_thread() != threading.main_thread()):
        # Silently ignore this error (expected when GUI is closed)
        return
    # For all other errors, show them normally
    original_excepthook(exc_type, exc_value, traceback)

sys.excepthook = suppress_thread_errors

# Now we can import and run the GUI

if __name__ == "__main__":
    try:
        from gui.main import main
        main()
    except Exception as e:
        print(f"Error launching GUI: {e}")
        import traceback
        traceback.print_exc()
        # Keep console open so user can see the error
        input("\nPress Enter to exit...")
        sys.exit(1)
