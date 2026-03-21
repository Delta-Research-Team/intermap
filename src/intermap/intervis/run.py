"""
Entry point for the InterMap Visualizations app.
Launches the application with the appropriate configuration.

"""

import atexit
import os
import signal
import socket
import sys
import threading
import time
import webbrowser
from pathlib import Path

import uvicorn
from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from shiny import App

ROOT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT_DIR))


STATIC_DIR = ROOT_DIR / "statics"
FAVICON_PATH = STATIC_DIR / "favicon-32x32.png"

_server_instance = None


def is_port_in_use(port):
    """Check if a port is already in use."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('127.0.0.1', port)) == 0


def kill_process_on_port(port):
    """Kill any process occupying the given port."""
    import platform

    system = platform.system()
    try:
        if system == "Windows":
            import subprocess
            result = subprocess.run(
                ['netstat', '-ano'], capture_output=True, text=True
            )
            for line in result.stdout.strip().split('\n'):
                if f'127.0.0.1:{port}' in line or f'0.0.0.0:{port}' in line:
                    parts = line.split()
                    pid = parts[-1]
                    if pid.isdigit() and int(pid) != os.getpid():
                        subprocess.run(['taskkill', '/F', '/PID', pid],
                                       capture_output=True)
                        print(f"Killed previous process on port {port} (PID: {pid})")
        else:
            import subprocess
            result = subprocess.run(
                ['lsof', '-ti', f':{port}'], capture_output=True, text=True
            )
            pids = result.stdout.strip().split('\n')
            for pid in pids:
                if pid.isdigit() and int(pid) != os.getpid():
                    os.kill(int(pid), signal.SIGKILL)
    except Exception as e:
        print(f"Warning: Could not kill process on port {port}: {e}")


def open_browser(port):
    """Opens the default web browser to the application URL."""
    webbrowser.open(f'http://127.0.0.1:{port}')


def create_app():
    """Creates and configures the Shiny application with FastAPI wrapper."""
    from app.main import app_ui, server
    fast_app = FastAPI()

    fast_app.mount("/static", StaticFiles(directory=str(STATIC_DIR)),
                   name="static")

    @fast_app.get('/favicon.ico')
    async def favicon():
        """Serve favicon.ico"""
        return FileResponse(str(FAVICON_PATH))

    @fast_app.get('/healthcheck')
    async def healthcheck():
        """Endpoint used by the browser to signal it's still alive."""
        return {"status": "alive"}
    app = App(ui=app_ui, server=server)
    fast_app.mount("/", app)

    return fast_app


class ServerWithShutdown(uvicorn.Server):
    """Custom uvicorn server that exposes a clean shutdown method."""

    def install_signal_handlers(self):
        """Override to handle signals ourselves."""
        signal.signal(signal.SIGINT, self._handle_signal)
        signal.signal(signal.SIGTERM, self._handle_signal)

    def _handle_signal(self, signum, frame):
        """Handle interrupt signals for clean shutdown."""
        self.should_exit = True


class BrowserMonitor:
    """
    Monitors browser connectivity via periodic HTTP checks.
    Shuts down the server when the browser is no longer reachable.
    """

    def __init__(self, port, server_instance, check_interval=5,
                 max_failures=3):
        self.port = port
        self.server = server_instance
        self.check_interval = check_interval
        self.max_failures = max_failures
        self._failures = 0
        self._running = True
        self._started = False
        self._thread = threading.Thread(target=self._monitor, daemon=True)

    def start(self):
        """Start monitoring after a grace period for initial startup."""
        self._thread.start()

    def stop(self):
        """Stop the monitor."""
        self._running = False

    def _monitor(self):
        """Periodically check if the browser is still connected."""
        import urllib.request
        import urllib.error

        time.sleep(10)
        self._started = True


        while self._running:
            time.sleep(self.check_interval)
            try:
                req = urllib.request.Request(
                    f'http://127.0.0.1:{self.port}/healthcheck',
                    method='GET'
                )
                urllib.request.urlopen(req, timeout=3)
                self._failures = 0
            except (urllib.error.URLError, ConnectionError, OSError):
                self._failures += 1
                if self._failures >= self.max_failures:
                    print(f"\nBrowser disconnected ({self._failures} failed "
                          f"checks). Shutting down server...")
                    self.server.should_exit = True
                    break


def shutdown_cleanup(port):
    """Cleanup function called at exit."""
    print(f"Cleaning up port {port}...")
    if is_port_in_use(port):
        kill_process_on_port(port)
    print("Shutdown complete.")


def main():
    """Main function to run the application."""
    port = 8000

    if is_port_in_use(port):
        kill_process_on_port(port)
        time.sleep(1)

        if is_port_in_use(port):
            print(f"ERROR: Could not free port {port}. Please close the "
                  f"process manually or use a different port.")
            sys.exit(1)
        print(f"Port {port} is now free.")

    atexit.register(shutdown_cleanup, port)
    threading.Timer(1.5, open_browser, args=[port]).start()

    print("Starting InterMap Visualizations...")
    print("Server Configuration:")
    print(f"- Host: 127.0.0.1")
    print(f"- Port: {port}")
    print(f"- Static Directory: {STATIC_DIR}")
    print(f"- Favicon Path: {FAVICON_PATH}")
    print(f"- Watch Directory: {ROOT_DIR}")
    print("- Auto-shutdown: Enabled (on browser disconnect)")

    config = uvicorn.Config(
        "run:create_app",
        host="127.0.0.1",
        port=port,
        reload=True,
        reload_dirs=[str(ROOT_DIR)],
        factory=True,
        log_level="info"
    )

    server = ServerWithShutdown(config)
    monitor = BrowserMonitor(port, server)
    monitor.start()
    server.run()
    monitor.stop()
    print("Server stopped. Goodbye!")


if __name__ == "__main__":
    main()
