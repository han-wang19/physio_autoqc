"""Start one experiment's manual review UI on its registered local port."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
import project_config as cfg


def review_command(experiment, port=None):
    selected_port = port if port is not None else cfg.EXPERIMENTS[experiment]['REVIEW_PORT']
    if not 1 <= selected_port <= 65535:
        raise ValueError('Port must be between 1 and 65535')
    command = [sys.executable, '-m', 'streamlit', 'run', str(PROJECT_ROOT / 'review_app.py'),
               '--server.address', '127.0.0.1', '--server.port', str(selected_port),
               '--server.headless', 'true', '--browser.gatherUsageStats', 'false']
    environment = os.environ.copy()
    environment['PHYSIO_EXPERIMENT'] = experiment
    environment.setdefault('NUMBA_CACHE_DIR', str(PROJECT_ROOT / '.cache' / 'numba'))
    return command, environment, selected_port


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', choices=tuple(cfg.EXPERIMENTS), required=True)
    parser.add_argument('--port', type=int)
    args = parser.parse_args()
    command, environment, port = review_command(args.experiment, args.port)
    print(f'{args.experiment}: http://127.0.0.1:{port}', flush=True)
    try:
        return subprocess.call(command, cwd=PROJECT_ROOT, env=environment)
    except KeyboardInterrupt:
        return 130


if __name__ == '__main__':
    raise SystemExit(main())
