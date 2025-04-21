"""
Entry point for the InterMap Visualizations app.
Launches the application with the appropriate configuration.

"""

import webbrowser
import threading
import uvicorn
import sys
from pathlib import Path
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from shiny import App

# Configure root directory and system path
ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))

# Configure static directories and paths
STATIC_DIR = ROOT_DIR / "statics"
FAVICON_PATH = STATIC_DIR / "image" / "favicon-32x32.png"

def open_browser(port):
    """Opens the default web browser to the application URL."""
    webbrowser.open(f'http://127.0.0.1:{port}')

def create_app():
    """Creates and configures the Shiny application with FastAPI wrapper."""
    from app.main import app_ui, server

    # Create FastAPI app
    fast_app = FastAPI()

    # Mount static files
    fast_app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    # Add favicon.ico endpoint
    @fast_app.get('/favicon.ico')
    async def favicon():
        """Serve favicon.ico"""
        return FileResponse(str(FAVICON_PATH))

    # Create and mount Shiny app
    app = App(ui=app_ui, server=server)
    fast_app.mount("/", app)

    return fast_app

def main():
    """Main function to run the application."""
    port = 8013

    # Start browser in a new thread after a short delay
    threading.Timer(1.5, open_browser, args=[port]).start()

    # Print configuration information
    print("=" * 50)
    print("Starting InterMap Visualizations...")
    print("-" * 50)
    print("Server Configuration:")
    print(f"- Host: 127.0.0.1")
    print(f"- Port: {port}")
    print(f"- Static Directory: {STATIC_DIR}")
    print(f"- Favicon Path: {FAVICON_PATH}")
    print(f"- Watch Directory: {ROOT_DIR}")
    print("- Reload Mode: Enabled")
    print("=" * 50)

    # Configure and run the application using uvicorn
    uvicorn.run(
        "run:create_app",
        host="127.0.0.1",
        port=port,
        reload=True,
        reload_dirs=[str(ROOT_DIR)],
        factory=True,
        log_level="info"
    )

if __name__ == "__main__":
    main()
