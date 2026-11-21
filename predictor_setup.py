"""Configure existing SignalP 4 Conda and MHCscan Phobius installations."""
from copy import deepcopy
import json
from pathlib import Path
import re
import shutil
import subprocess


def configure(args, config):
    result = deepcopy(config)
    if args.signalp != 'skip':
        settings = result.setdefault('signalp', {})
        settings.setdefault('format', 'signalp'+args.signalp_version)
        if args.signalp4_env:
            if args.signalp_version != '4': raise ValueError('--signalp4-env requires --signalp-version 4')
            if 'argv' in settings: raise ValueError('Use --signalp4-env or config signalp argv, not both')
            conda = shutil.which(args.conda_executable)
            if not conda: raise ValueError(f'{args.conda_executable} not found; use --conda-executable /home/ark/miniconda3/bin/conda')
            env = args.signalp4_env
            if Path(env).is_dir(): prefix = Path(env).resolve()
            else:
                probe = subprocess.run([conda, 'env', 'list', '--json'], capture_output=True, text=True, timeout=60)
                if probe.returncode: raise ValueError(f'Cannot list Conda environments: {probe.stderr.strip()}')
                try: matches = [Path(p) for p in json.loads(probe.stdout)['envs'] if Path(p).name==env]
                except (ValueError, KeyError) as exc: raise ValueError('Conda env list did not return valid JSON') from exc
                if len(matches)!=1: raise ValueError(f'Conda environment {env!r} is missing or ambiguous; use its full prefix with --signalp4-env')
                prefix = matches[0].resolve()
            script = prefix/'bin/signalp'
            if not script.is_file(): raise ValueError(f'{script} not found in requested SignalP 4 environment')
            launcher = [conda, 'run', '--no-capture-output', '-p', str(prefix), 'signalp']
            settings.update(argv=launcher+['-t', '{signalp_gram}', '-f', 'short', '-s', 'best', '-c', '70', '{input}'],
                            output='stdout', version_argv=launcher+['-V'], format='signalp4')
            settings['model_files'] = list(dict.fromkeys(settings.get('model_files', [])+[str(script.resolve())]))
            # The script checksum is known; hidden network/model files remain user-configurable.
            settings['model_provenance_note'] = 'SignalP executable tracked; add installed network/model files to config for complete model hashes'
    if args.phobius != 'skip' and args.phobius_dir:
        settings = result.setdefault('phobius', {})
        if 'argv' in settings: raise ValueError('Use --phobius-dir or config phobius argv, not both')
        folder = args.phobius_dir.resolve()
        paths = [folder/name for name in ['phobius.pl', 'decodeanhmm.64bit', 'phobius.options', 'phobius.model']]
        for path in paths:
            if not path.is_file(): raise ValueError(f'Incomplete MHCscan Phobius installation: missing {path}')
        if not paths[1].stat().st_mode & 0o111:
            raise ValueError(f'Phobius decoder is not executable. Run: chmod +x {paths[1]}')
        perl = shutil.which('perl')
        if not perl: raise ValueError('Perl is required for the MHCscan Phobius script')
        # Its internal shell command does not quote the decoder/model directory.
        if re.search(r'[^a-zA-Z0-9_./-]', str(folder)):
            raise ValueError('MHCscan Phobius requires an installation directory without spaces or shell metacharacters; move the complete bundle to a simple path')
        settings.update(argv=[perl, str(paths[0]), '-short'], output='stdout', stdin=True)
        settings['model_files'] = list(dict.fromkeys(settings.get('model_files', [])+[str(p) for p in paths]))
        match = re.search(r'phobius\.pl\s+([\d.]+)', paths[0].read_text(errors='replace')[:500])
        if match: settings.setdefault('version', 'Phobius '+match[1])
    return result
