"""Read-only audit of the files Git would include in a public code snapshot."""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GIT = ['git', '-c', 'safe.directory='+ROOT.as_posix()]
DATA_SUFFIXES = {'.fif','.bdf','.edf','.npz','.npy','.csv','.tsv','.xlsx','.pptx','.png','.jpg','.jpeg',
                 '.mat','.h5','.hdf5','.set','.fdt','.pkl','.pickle','.ipynb'}
EXAMPLES = {'examples/events.csv', 'examples/events.tsv'}
PRIVATE_DIRS = {'outputs','data','.cache','.local_docs','.venv','venv','__pycache__'}
PRIVATE_FILES = {'project_config.local.json','project_config.local.py','try.py'}
DRIVE_PATH = re.compile(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]')


def public_candidates():
    result = subprocess.run(GIT+['ls-files','--cached','--others','--exclude-standard','-z'],
                            cwd=ROOT, check=True, stdout=subprocess.PIPE)
    return sorted(set(name.decode('utf-8') for name in result.stdout.split(b'\0') if name))


def audit_files(root, names):
    problems = []
    for name in names:
        path = root/name
        if not path.is_file():
            continue  # Tracked deletions will not be present in the next snapshot.
        parts = Path(name).parts
        if (any(part in PRIVATE_DIRS for part in parts) or Path(name).name in PRIVATE_FILES
                or (Path(name).name.startswith('.env') and Path(name).name != '.env.example')):
            problems.append((name,'private file is included by Git'))
            continue
        if (path.suffix.lower() in DATA_SUFFIXES or name.lower().endswith('.fif.gz')) and name not in EXAMPLES:
            problems.append((name,'data/result file is included by Git'))
            continue
        if path.suffix.lower() in {'.py','.md','.json','.txt','.yml','.yaml'}:
            if DRIVE_PATH.search(path.read_text(encoding='utf-8', errors='replace')):
                problems.append((name,'machine-specific absolute Windows path'))
    return problems


def main():
    names = public_candidates()
    problems = audit_files(ROOT,names)
    for name,reason in problems:
        print(f'FAIL: {name}: {reason}')
    if problems:
        print('Public snapshot check failed. No files were changed or uploaded.')
        return 1
    print(f'Public snapshot check passed: {sum((ROOT/n).is_file() for n in names)} code/document files.')
    print('No study data, private configuration or machine-specific Windows paths found.')
    print('This check covers the current candidate files, not historical commit contents.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
