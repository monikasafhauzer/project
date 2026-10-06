import os
import shutil
import subprocess


def dependencies():
    result = {}
    for tool in ('minimap2', 'vsearch'):
        path = shutil.which(tool)
        if not path:
            result[tool] = 'Missing — install and add to PATH'
            continue
        try:
            proc = subprocess.run([path, '--version'], capture_output=True, text=True, timeout=10)
            result[tool] = (proc.stdout + proc.stderr).strip().splitlines()[0] if proc.returncode == 0 else 'Unable to execute'
        except (OSError, subprocess.TimeoutExpired):
            result[tool] = 'Unable to execute'
    return result


def run(command, timeout=900):
    if not shutil.which(command[0]):
        raise RuntimeError(f'{command[0]} is missing. Install it and add its executable to PATH, then restart the app.')
    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'{command[0]} exceeded {timeout} seconds. Try fewer reads or run on a larger machine.') from exc
    if proc.returncode:
        raise RuntimeError(f'{command[0]} failed: {proc.stderr[-4000:]}')
    return proc.stdout


def threads():
    return str(min(4, os.cpu_count() or 1))
