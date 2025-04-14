"""
Entry point for the InterMap Visualizations app.
Launches the application with the appropriate configuration.

"""

import sys
import threading
import webbrowser
from pathlib import Path

import uvicorn

ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))


def open_browser(port):
    """Opens the default web browser to the application URL."""
    webbrowser.open(f'http://127.0.0.1:{port}')


def main():
    """Main function to run the application."""
    port = 8000

    # Start browser in a new thread after a short delay
    threading.Timer(1.5, open_browser, args=[port]).start()

    # Configure and run the application using uvicorn
    uvicorn.run(
        "app.main:app", host="127.0.0.1", port=port, reload=True,
        reload_dirs=[str(ROOT_DIR)],
    )


if __name__ == "__main__":
    main()
