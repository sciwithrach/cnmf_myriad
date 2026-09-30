"""Record how an analysis was run, next to its outputs: analysis_call.txt (to read) and analysis_call.json (to parse)."""
from __future__ import annotations

import json
import platform
import shlex
import time
from importlib import metadata
from pathlib import Path
from typing import Dict, Optional

PACKAGES = ('scanpy', 'anndata', 'decoupler', 'pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn')


def resolved_command(script: str, settings: Dict[str, object]) -> str:
    """The command line with every option spelled out, defaults included (settings: option name -> value)."""
    parts = [script]
    for key, value in settings.items():
        if value is None or value is False:
            continue
        parts.append(f'--{key}' if value is True else f'--{key} {shlex.quote(str(value))}')
    return ' '.join(parts)


def describe_input(path: str | Path) -> Dict[str, object]:
    """Absolute path (symlinks resolved), size and modification time of an input file or folder."""
    p = Path(path).resolve()
    if not p.exists():
        return {'path': str(p), 'exists': False}
    stat = p.stat()
    return {'path': str(p), 'exists': True, 'kind': 'folder' if p.is_dir() else 'file',
            'size_bytes': None if p.is_dir() else stat.st_size,
            'modified': time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(stat.st_mtime))}


def package_versions() -> Dict[str, str]:
    """Python and the versions of the packages the analysis uses."""
    versions = {'python': platform.python_version()}
    for name in PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = 'not installed'
    return versions


def write_call_record(out: str | Path, script: str, settings: Dict[str, object], inputs: Dict[str, str],
                      typed: Optional[str] = None, nextflow_info: Optional[str | Path] = None) -> None:
    """Write analysis_call.txt and analysis_call.json into `out`.

    settings : every option with its final value (defaults included)
    inputs : name -> path of the input files and folders
    typed : the command line as it was typed, if different from the resolved one
    nextflow_info : optional text file (one 'key: value' per line) describing the Nextflow run that started this
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    when = time.strftime('%Y-%m-%d %H:%M:%S')
    command = resolved_command(script, settings)
    described = {name: describe_input(path) for name, path in inputs.items()}
    versions = package_versions()

    nextflow: Dict[str, str] = {}
    if nextflow_info and Path(nextflow_info).exists():
        for line in Path(nextflow_info).read_text().splitlines():
            key, _, value = line.partition(':')
            if value:
                nextflow[key.strip()] = value.strip()

    record = {'date': when, 'command': command, 'command_as_typed': typed, 'settings': settings,
              'inputs': described, 'versions': versions, 'nextflow': nextflow}
    (out / 'analysis_call.json').write_text(json.dumps(record, indent=2, default=str) + '\n')

    lines = ['GEP analysis call', '=================', f'Date: {when}', '']
    if nextflow:
        lines += ['Nextflow run'] + [f'  {k}: {v}' for k, v in nextflow.items()] + ['']
    lines += ['Command, with every option and default written out', f'  {command}', '']
    if typed and typed != command:
        lines += ['Command as typed', f'  {typed}', '']
    lines += ['Settings'] + [f'  {k}: {v}' for k, v in settings.items()] + ['']
    lines += ['Inputs']
    for name, info in described.items():
        detail = 'not found' if not info['exists'] else (
            f"{info['kind']}, modified {info['modified']}" + (f", {info['size_bytes']:,} bytes" if info['size_bytes'] else ''))
        lines += [f'  {name}: {info["path"]}  ({detail})']
    lines += ['', 'Versions'] + [f'  {k}: {v}' for k, v in versions.items()]
    (out / 'analysis_call.txt').write_text('\n'.join(lines) + '\n')
