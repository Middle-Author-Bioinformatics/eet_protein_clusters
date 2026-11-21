"""Native predictor adapters and conservative, source-preserving refinements.

Only local processes or explicitly supplied files are used. No web submissions,
licensed downloads, or inference of extracellular exposure from an SP alone.
"""
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from common import LOG, mark_stage, sha256, signature, stage_reusable

TOOLS = ('psortb', 'signalp', 'phobius', 'lipop')
EXTENSIONS = dict(psortb='psortb.tsv', signalp='signalp.tsv', phobius='phobius.txt', lipop='lipop.txt')
EXECUTABLES = dict(psortb='psort', signalp='signalp6', phobius='phobius.pl', lipop='LipoP')


def parse_native(tool, path):
    """Map native first-token IDs to feature dictionaries; reject duplicates/invalid data."""
    result = {}
    def add(pid, row):
        pid = pid.strip().split()[0].lstrip('>')
        if pid in result: raise ValueError(f'{path}: duplicate {tool} prediction for {pid}')
        result[pid] = row
    lines = Path(path).read_text().splitlines()
    signalp4 = tool=='signalp' and any(re.match(r'#\s+(?:SignalP-4|name\s+Cmax)', line) for line in lines)
    if signalp4:
        for line in lines:
            if not line.strip() or line.lstrip().startswith('#'): continue
            words = line.split()
            if len(words)!=12 or words[9] not in ['Y', 'N'] or words[11] not in ['SignalP-TM', 'SignalP-noTM']:
                raise ValueError(f'{path}: expected native SignalP 4 -f short output (12 columns)')
            try:
                values = [float(words[i]) for i in [1, 3, 5, 7, 8, 10]]
                positions = [int(words[i]) for i in [2, 4, 6]]
                if not all(math.isfinite(x) and 0<=x<=1 for x in values) or min(positions)<1: raise ValueError()
            except ValueError as exc: raise ValueError(f'{path}: invalid SignalP 4 scores/positions') from exc
            # Ymax is the FIRST mature residue, not the last SP residue.
            mature_start = positions[1]
            if words[9]=='Y' and mature_start<2: raise ValueError(f'{path}: invalid positive SignalP 4 cleavage position')
            add(words[0], dict(prediction='SP' if words[9]=='Y' else 'OTHER', generation='4',
                               cleavage=f'{mature_start-1}-{mature_start}' if words[9]=='Y' else '',
                               d_score=words[8], d_cutoff=words[10], networks=words[11],
                               lipoprotein_status='not_assessed_by_SignalP4', native=dict(zip(
                                   ['id','Cmax','Cmax_pos','Ymax','Ymax_pos','Smax','Smax_pos','Smean','D','SP','Dcutoff','networks'], words))))
    elif tool in ['psortb', 'signalp']:
        key = 'SeqID' if tool == 'psortb' else 'ID'
        # Accept spelling used by older PSORTb releases too.
        header_index = next((i for i, line in enumerate(lines)
                             if line.lstrip('# ').split('\t')[0].lower() == key.lower()), None)
        if header_index is None: raise ValueError(f'{path}: missing {tool} native TSV header ({key})')
        rows = [lines[header_index].lstrip('# ')] + [s for s in lines[header_index+1:] if s.strip() and not s.startswith('#')]
        reader = csv.DictReader(rows, delimiter='\t')
        for row in reader:
            row = {k.strip(): (v.strip() if v is not None else '') for k, v in row.items() if k is not None}
            id_key = next((k for k in row if k.lower() == key.lower()), None)
            if tool == 'psortb':
                loc = row.get('Final_Localization', row.get('Localization', ''))
                if not loc: raise ValueError(f'{path}: PSORTb needs Final_Localization (long) or Localization (terse)')
                add(row[id_key], dict(localization=loc, score=row.get('Final_Score', row.get('Score', '')),
                                      secondary=row.get('Secondary_Localization', ''), native=row))
            else:
                prediction = row.get('Prediction', '')
                if prediction not in ['OTHER', 'SP', 'LIPO', 'TAT', 'TATLIPO', 'PILIN']:
                    raise ValueError(f'{path}: unsupported SignalP 6 prediction {prediction!r}')
                add(row[id_key], dict(prediction=prediction, generation='6', cleavage=row.get('CS Position', ''), lipoprotein_status='assessed', native=row))
    elif tool == 'phobius':
        for line in lines:
            words = line.split()
            if not words or words[0].startswith('#') or words[0] in ['SEQENCE', 'SEQUENCE', 'ID']: continue
            if len(words) != 4 or not words[1].isdigit() or words[2] not in ['Y', '0']:
                raise ValueError(f'{path}: expected Phobius -short output: ID TM SP PREDICTION')
            add(words[0], dict(tm_count=int(words[1]), signal_peptide=words[2]=='Y', topology=words[3]))
    elif tool == 'lipop':
        for line in lines:
            m = re.match(r'^#\s+(\S+)\s+(SpII|SpI|TMH|CYT)\s+score=(\S+)\s+margin=(\S+)(.*)$', line)
            if m:
                add(m[1], dict(prediction=m[2], score=m[3], margin=m[4], detail=m[5].strip()))
    if not result: raise ValueError(f'{path}: no recognizable {tool} predictions')
    return result


def map_predictions(predictions, mapping):
    """Prefer safe IDs; accept unambiguous original IDs only within this genome."""
    allowed = {r['protein_id'] for r in mapping}
    originals = {}
    for row in mapping: originals.setdefault(row['original_id'], []).append(row['protein_id'])
    result = {}
    for pid, values in predictions.items():
        if pid in allowed: target = pid
        elif pid in originals and len(originals[pid]) == 1: target = originals[pid][0]
        else: raise ValueError(f'Unknown or ambiguous localization protein ID {pid!r}; use IDs from inputs/proteomes/*.faa')
        if target in result: raise ValueError(f'Conflicting localization IDs map to {target}')
        result[target] = values
    return result


def default_command(tool, input_path, outdir, gram, threads, signalp_mode, signalp_format='signalp6'):
    if tool == 'psortb': return ['psort', '-p' if gram=='positive' else '-n', '-o', 'long', str(input_path)], 'stdout'
    if tool == 'signalp' and signalp_format=='signalp4':
        return ['signalp', '-t', 'gram+' if gram=='positive' else 'gram-', '-f', 'short', '-s', 'best', '-c', '70', str(input_path)], 'stdout'
    if tool == 'signalp':
        return ['signalp6', '--fastafile', str(input_path), '--organism', 'other', '--output_dir', str(outdir),
                '--format', 'none', '--mode', signalp_mode, '--torch_num_threads', str(threads), '--write_procs', '1'], 'prediction_results.txt'
    if tool == 'phobius': return ['phobius.pl', '-short'], 'stdout'
    raise ValueError('LipoP execution requires its installation-specific argv in --localization-config; native short output can also be imported')


def collect(genome, mapping, args, config, output, manifest):
    """Run/import selected analyses with per-genome checksum caches and exact ID validation."""
    data = {}; audit = []
    for tool in TOOLS:
        settings = config.get(tool, {})
        mode = getattr(args, tool)
        imported = args.localization_import / f'{genome["genome"]}.{EXTENSIONS[tool]}' if args.localization_import else None
        dest = output/f'{tool}.txt'; folder = output/tool
        entry = dict(genome=genome['genome'], tool=tool, status='unavailable', reason='', proteins_assessed=0)
        if mode == 'skip': entry['reason'] = 'explicitly skipped'
        elif imported and imported.exists():
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(imported, dest)
            entry.update(status='imported', input=str(imported.resolve()), sha256=sha256(dest),
                         version=settings.get('version', 'unavailable (not supplied)'),
                         model_files={str(Path(p).resolve()): sha256(p) for p in settings.get('model_files', [])})
        else:
            if tool == 'lipop' and mode == 'auto' and 'argv' not in settings:
                entry['reason'] = 'optional LipoP not configured; SignalP 6 provides lipoprotein predictions when available'
                audit.append(entry); data[tool] = {}; continue
            defaults = None
            if 'argv' not in settings:
                selected_executable = 'signalp' if tool=='signalp' and settings.get('format')=='signalp4' else EXECUTABLES[tool]
                if mode == 'auto' and not shutil.which(selected_executable):
                    entry['reason'] = f'{selected_executable} not found; no import supplied'
                    audit.append(entry); data[tool] = {}; continue
                defaults = default_command(tool, genome['fasta'], folder.resolve(), genome['gram'], args.threads, args.signalp_mode, settings.get('format', 'signalp6'))
            template = settings.get('argv', defaults[0] if defaults else [])
            if not isinstance(template, list) or not template or not all(isinstance(s, str) for s in template):
                raise ValueError(f'{tool}: config argv must be a nonempty list of strings (no shell command)')
            context = dict(input=str(Path(genome['fasta']).resolve()), outdir=str(folder.resolve()),
                           signalp_gram='gram+' if genome['gram']=='positive' else 'gram-',
                           gram=genome['gram'], gram_flag='-p' if genome['gram']=='positive' else '-n', threads=args.threads)
            try: argv = [s.format(**context) for s in template]
            except KeyError as exc: raise ValueError(f'{tool}: unknown argv placeholder {exc}') from exc
            binary = shutil.which(argv[0])
            if not binary: raise ValueError(f'{tool}: {argv[0]} not found. Install/configure the predictor or provide --localization-import; use --{tool} skip to omit it explicitly. See GENOME_WORKFLOW.md.')
            version = settings.get('version')
            if not version and tool in ['psortb', 'signalp']:
                probe_argv = settings.get('version_argv', [binary, '-V' if tool=='signalp' and settings.get('format')=='signalp4' else '--version'])
                manifest['commands'].append(probe_argv)
                probe = subprocess.run(probe_argv, capture_output=True, text=True, timeout=30)
                version = (probe.stdout+'\n'+probe.stderr).strip() if probe.returncode == 0 else None
            version = version or 'unavailable (supply version in config)'
            model_hashes = {str(Path(p).resolve()): sha256(p) for p in settings.get('model_files', [])}
            use_stdin = settings.get('stdin', tool=='phobius' and 'argv' not in settings)
            result_name = settings.get('output', defaults[1] if defaults else 'stdout')
            provenance = dict(argv=argv, stdin=use_stdin, result_file=result_name, result_format=settings.get('format', tool), executable=str(Path(binary).resolve()), executable_sha256=sha256(binary),
                              version=version, model_files=model_hashes, cache_tag=args.predictor_cache_tag)
            if settings.get('model_provenance_note'): manifest['warnings'].append(settings['model_provenance_note'])
            if not model_hashes:
                manifest['warnings'].append(f'{genome["genome"]}/{tool}: model checksums unavailable; supply config model_files or change --predictor-cache-tag after updating models.')
            key = signature(dict(input=sha256(genome['fasta']), **provenance))
            stamp = output/f'{tool}.stage.json'
            output.mkdir(parents=True, exist_ok=True); folder.mkdir(exist_ok=True)
            reused = stage_reusable(stamp, key, [dest], args.resume and not args.force)
            if not reused:
                result_path = folder/result_name
                if result_name != 'stdout':
                    if Path(result_name).is_absolute() or '..' in Path(result_name).parts:
                        raise ValueError(f'{tool}: output must be a relative filename beneath {{outdir}} or stdout')
                    result_path.unlink(missing_ok=True)
                logfile = output/f'{tool}.run.log'
                run_prediction(argv, genome['fasta'] if use_stdin else None, dest if result_name=='stdout' else folder/'stdout.log', logfile, manifest)
                if result_name == 'stdout': pass
                elif result_path.is_file(): shutil.copyfile(result_path, dest)
                else: raise ValueError(f'{tool}: command succeeded but expected output {result_path} is missing; check localization config')
            entry.update(status='ran', cache_reused=reused, sha256=sha256(dest), **provenance)
        if entry['status'] != 'unavailable':
            data[tool] = map_predictions(parse_native(tool, dest), mapping)
            if tool=='signalp': entry['detected_generations'] = sorted({v.get('generation', 'unknown') for v in data[tool].values()})
            entry['proteins_assessed'] = len(data[tool])
            if len(data[tool]) < len(mapping):
                manifest['warnings'].append(f'{genome["genome"]}/{tool}: {len(mapping)-len(data[tool])} proteins have no result; labeled unavailable.')
            if entry['status'] == 'ran': mark_stage(output/f'{tool}.stage.json', key, [dest])
        else: data[tool] = {}
        audit.append(entry)
    manifest['localization'].extend(audit)
    manifest['skipped_analyses'].extend(e for e in audit if e['status']=='unavailable')
    return data


def run_prediction(argv, input_path, stdout_path, logfile, manifest):
    """Keep stderr banners/warnings out of machine-readable stdout results."""
    from contextlib import ExitStack
    LOG.info('Running %s', argv[0]); manifest['commands'].append(argv)
    with ExitStack() as stack:
        source = stack.enter_context(open(input_path, 'rb')) if input_path else subprocess.DEVNULL
        out = stack.enter_context(open(stdout_path, 'w'))
        err = stack.enter_context(open(logfile, 'w'))
        proc = subprocess.run(argv, stdin=source, stdout=out, stderr=err)
    if proc.returncode:
        tail = '\n'.join(Path(logfile).read_text(errors='replace').splitlines()[-12:])
        raise ValueError(f'{argv[0]} exited {proc.returncode}. See {logfile}\n{tail}')


LOCALIZATION_FIELDS = ['psortb_status', 'psortb_localization', 'psortb_score', 'psortb_secondary',
                       'signalp_status', 'signalp_prediction', 'signalp_generation', 'signalp_cleavage',
                       'signalp_d_score', 'signalp_d_cutoff', 'signalp_networks', 'lipoprotein_assessment', 'phobius_status',
                       'phobius_tm_count', 'phobius_signal_peptide', 'phobius_topology', 'lipop_status',
                       'lipop_prediction', 'secretion_evidence', 'membrane_attachment_evidence',
                       'cell_wall_evidence', 'extracellular_exposure_evidence', 'localization_interpretation',
                       'inference_basis', 'localization_conflict']


def infer(evidence, gram):
    """Refine Unknown only to feature-supported categories; never silently override PSORTb."""
    p, s, t, l = (evidence.get(tool, {}) for tool in TOOLS)
    baseline = p.get('localization', 'unavailable')
    sp = s.get('prediction')
    secretion = []
    if sp and sp != 'OTHER': secretion.append(f'SignalP:{sp}')
    if t.get('signal_peptide'): secretion.append('Phobius:SP')
    if l.get('prediction') in ['SpI', 'SpII']: secretion.append('LipoP:'+l['prediction'])
    lipid = sp in ['LIPO', 'TATLIPO'] or l.get('prediction') == 'SpII'
    tm = t.get('tm_count')  # absence is not zero
    membrane = []
    if lipid: membrane.append('predicted_lipid_anchor')
    if tm is not None and tm > 0: membrane.append(f'Phobius:{tm}_TM')
    if baseline in ['CytoplasmicMembrane', 'OuterMembrane']: membrane.append('PSORTb:'+baseline)
    conflicts = []
    if baseline == 'Cytoplasmic' and (secretion or membrane): conflicts.append('PSORTb_cytoplasmic_vs_export_or_attachment')
    if sp == 'OTHER' and (t.get('signal_peptide') or l.get('prediction') in ['SpI', 'SpII']): conflicts.append('signal_peptide_predictors_disagree')
    if sp and sp != 'OTHER' and t and not t.get('signal_peptide'): conflicts.append('signal_peptide_predictors_disagree')
    if baseline == 'Extracellular' and membrane: conflicts.append('PSORTb_extracellular_vs_attachment')
    inferred = baseline; basis = 'PSORTb' if p else 'no_PSortb_result'
    if baseline.lower() in ['unknown', 'unavailable']:
        if lipid: inferred = 'putative_membrane_attached_lipoprotein'; basis = 'SignalP_and/or_LipoP'
        elif tm is not None and tm > 0: inferred = 'putative_membrane_protein'; basis = 'Phobius_TM'
        elif secretion:
            inferred = 'export_signal_detected_destination_unresolved'; basis = 'signal_peptide_prediction'
        else:
            inferred = 'Unknown' if p else ('unresolved_no_positive_targeting_evidence' if s or t or l else 'unavailable')
            basis = 'insufficient_evidence'
    return dict(psortb_status='assessed' if p else 'unavailable', psortb_localization=baseline,
                psortb_score=p.get('score', ''), psortb_secondary=p.get('secondary', ''),
                signalp_status='assessed' if s else 'unavailable', signalp_prediction=sp or 'unavailable',
                signalp_generation=s.get('generation', 'unavailable'), signalp_cleavage=s.get('cleavage', ''),
                signalp_d_score=s.get('d_score', ''), signalp_d_cutoff=s.get('d_cutoff', ''),
                signalp_networks=s.get('networks', ''),
                lipoprotein_assessment='assessed' if (sp in ['LIPO', 'TATLIPO'] or s.get('generation')=='6' or l) else 'unavailable', phobius_status='assessed' if t else 'unavailable',
                phobius_tm_count=tm if tm is not None else 'unavailable',
                phobius_signal_peptide=str(t['signal_peptide']).lower() if t else 'unavailable',
                phobius_topology=t.get('topology', ''), lipop_status='assessed' if l else 'unavailable',
                lipop_prediction=l.get('prediction', 'unavailable'),
                secretion_evidence=';'.join(secretion) or ('no_signal_predicted' if s or t or l else 'unavailable'),
                membrane_attachment_evidence=';'.join(membrane) or ('no_attachment_predicted' if t else 'unavailable'),
                cell_wall_evidence='PSORTb:'+baseline if baseline in ['Cellwall', 'CellWall'] else 'unresolved',
                extracellular_exposure_evidence='PSORTb:Extracellular' if baseline=='Extracellular' else
                    ('PSORTb:OuterMembrane_surface_exposure_unresolved' if baseline=='OuterMembrane' else 'unresolved'),
                localization_interpretation=inferred, inference_basis=basis,
                localization_conflict=';'.join(sorted(set(conflicts))))
