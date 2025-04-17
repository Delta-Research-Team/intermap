"""
Entry point for the InterMap Visualizations app.
Launches the application with the appropriate configuration.

"""

import webbrowser
import threading
import uvicorn

import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))

def open_browser(port):
    """Opens the default web browser to the application URL."""
    webbrowser.open(f'http://127.0.0.1:{port}')

def main():
    """Main function to run the application."""
    port = 8001

    # Start browser in a new thread after a short delay
    threading.Timer(1.5, open_browser, args=[port]).start()

    print("=" * 50)
    print("Starting InterMap Visualizations...")
    print("-" * 50)
    print("Server Configuration:")
    print(f"- Host: 127.0.0.1")
    print(f"- Port: {port}")
    print("- Reload Mode: Enabled")
    print(f"- Watch Directory: {ROOT_DIR}")
    print("=" * 50)

    # Configure and run the application using uvicorn
    uvicorn.run(
        "app.main:app",    # Import path to the application
        host="127.0.0.1",  # Listen only on localhost
        port=port,         # Use fixed port
        reload=True,       # Enable auto-reload on code changes
        reload_dirs=[str(ROOT_DIR)],  # Directories to watch for changes
    )

if __name__ == "__main__":
    main()
