#!/usr/bin/env python3
"""
Resubmit the jobs that have no output in hists/<processor>/ as one condor cluster.
See "Resubmit with condor" in README.md for all the options.

./resubmit_condor.py -p hadmonotop2022_0802 -m 2022_private_v2 -c lpc
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from optparse import OptionParser

# ---------------------------------------------------------------------------
# One JDL blob per cluster.  {sample} is the condor macro fed by the job list,
# every other {...} is filled in by this script.
# ---------------------------------------------------------------------------
JDL = {}

JDL['kisti'] = """universe                = vanilla
executable              = {executable}
arguments               = {arguments}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = {input_files}
{proxy_lines}output                  = logs/condor/run/out/{processor}_{sample}_$(Cluster)_$(Process).stdout
error                   = logs/condor/run/err/{processor}_{sample}_$(Cluster)_$(Process).stderr
log                     = logs/condor/run/log/{processor}_{sample}_$(Cluster)_$(Process).log
transfer_output_remaps  = "{processor}_{sample}.futures={outdir}/{sample}.futures"
initialdir              = {pwd}
JobBatchName            = {batchname}
accounting_group        = {accounting_group}
request_cpus            = {cpus}
request_memory          = {memory}
{extra_lines}queue SAMPLE from {listfile}
"""

JDL['lpc'] = """universe                = vanilla
executable              = {executable}
arguments               = {arguments}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = {input_files}
{proxy_lines}output                  = logs/condor/run/out/{processor}_{sample}_$(Cluster)_$(Process).stdout
error                   = logs/condor/run/err/{processor}_{sample}_$(Cluster)_$(Process).stderr
log                     = logs/condor/run/log/{processor}_{sample}_$(Cluster)_$(Process).log
transfer_output_remaps  = "{processor}_{sample}.futures={outdir}/{sample}.futures"
initialdir              = {pwd}
JobBatchName            = {batchname}
request_cpus            = {cpus}
request_memory          = {memory}
+SingularityImage       = "{image}"
{extra_lines}queue SAMPLE from {listfile}
"""

JDL['knut3'] = """universe                = vanilla
executable              = {executable}
arguments               = {arguments}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_input_files    = {input_files}
{proxy_lines}output                  = logs/condor/run/out/{processor}_{sample}_$(Cluster)_$(Process).stdout
error                   = logs/condor/run/err/{processor}_{sample}_$(Cluster)_$(Process).stderr
log                     = logs/condor/run/log/{processor}_{sample}_$(Cluster)_$(Process).log
transfer_output_remaps  = "{processor}_{sample}.futures={outdir}/{sample}.futures"
initialdir              = {pwd}
JobBatchName            = {batchname}
run_as_owner            = True
request_cpus            = {cpus}
request_memory          = {memory} MB
{extra_lines}queue SAMPLE from {listfile}
"""

# per cluster defaults, overridable with the options below or --site-config
SITE_DEFAULTS = {
    'kisti': {'cpus': '8', 'memory': '7000', 'accounting_group': 'group_cms', 'image': ''},
    'lpc':   {'cpus': '8', 'memory': '16000', 'accounting_group': '',
              'image': '/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base-almalinux8:0.7.22-py3.8'},
    'knut3': {'cpus': '8', 'memory': '2048', 'accounting_group': '', 'image': ''},
}

# Where each site keeps the tarballs that run.sh fetches with xrdcp.  {user} is
# filled in from the environment, so nothing here is tied to one account.
STORE = {
    'kisti': 'root://cms-xrdr.private.lo:2094//xrd/store/user/{user}',
    'lpc':   'root://cmseos.fnal.gov//store/user/{user}',
    'knut3': 'root://cluster142.knu.ac.kr//store/user/{user}',
}
# sites where xrdcp -f cannot overwrite, so the old file has to go first
XRDFS_RM = {'kisti': 'root://cms-xrdr.private.lo:2094/'}

# this script lives in <decaf>/analysis/, so the working copy is one level up;
# the tarballs are written next to it
DECAF_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARBALL_DIR = os.path.dirname(DECAF_DIR)

# tarball names the job script unpacks; keep in sync with it if you rename them
DECAF_TGZ = 'decaf.tgz'
PYLOCAL_TGZ = 'pylocal_3_8.tgz'

# keys that --site-config may set; same names as the long options below
SITE_KEYS = ('cluster', 'executable', 'cpus', 'memory', 'disk', 'image',
             'accounting_group', 'proxy', 'input_files', 'exec_args', 'extra',
             'retries', 'chunk', 'submit_opts')


def whoami():
    return os.environ.get('USER', os.environ.get('LOGNAME', 'unknown'))


def default_pylocal_dir():
    """site-packages directory of the environment, taken from $PYTHONPATH."""
    for path in os.environ.get('PYTHONPATH', '').split(':'):
        path = path.rstrip('/')
        if path.endswith('site-packages') and os.path.isdir(path):
            return path
    import site
    try:
        usersite = site.getusersitepackages()
    except AttributeError:
        usersite = ''
    if usersite and os.path.isdir(usersite):
        return usersite
    import sysconfig
    purelib = sysconfig.get_paths().get('purelib', '')
    return purelib if os.path.isdir(purelib) else ''


DRY_RUN = False          # set from -n/--dry-run, makes run() print only


def run(cmd, ignore_errors=False):
    """Run a command, printing it first; stop the script if it fails."""
    print('+ ' + ' '.join(cmd))
    if DRY_RUN:
        return
    rc = subprocess.call(cmd)
    if rc != 0 and not ignore_errors:
        sys.exit('command failed with exit code %d' % rc)


def tarball_paths():
    """(decaf tarball, pylocal tarball) that --tar writes and --copy sends."""
    return (os.path.join(TARBALL_DIR, DECAF_TGZ), os.path.join(TARBALL_DIR, PYLOCAL_TGZ))


def make_tarballs():
    """tar the decaf working copy and the python site-packages."""
    decaf_dir = DECAF_DIR
    decaf_tgz, pylocal_tgz = tarball_paths()

    excludes = ['analysis/logs', 'analysis/plots', 'analysis/datacards',
                'analysis/results', 'analysis/data/models',
                'analysis/hists/*/*.futures', 'analysis/hists/*/*.merged',
                'analysis/hists/*/*.reduced']
    cmd = ['tar', '--exclude-caches-all', '--exclude-vcs']
    cmd += ['--exclude=' + e for e in excludes]
    cmd += ['-czf', decaf_tgz, '-C', os.path.dirname(decaf_dir), os.path.basename(decaf_dir)]
    run(cmd)
    print('wrote ' + decaf_tgz)

    pydir = default_pylocal_dir()
    if not pydir:
        sys.exit('site-packages directory not found: check $PYTHONPATH')
    run(['tar', '--exclude-caches-all', '--exclude-vcs', '-czf', pylocal_tgz,
         '-C', os.path.dirname(pydir), os.path.basename(pydir)])
    print('wrote ' + pylocal_tgz)


def remote_path(url):
    """'root://host:port//store/user/me/x.tgz' -> '/store/user/me/x.tgz'"""
    match = re.match(r'^[a-z]+://[^/]+/+(/.*)$', url)
    if match:
        return match.group(1)
    return re.sub(r'^[a-z]+://[^/]+', '', url)


def copy_tarballs(options):
    """xrdcp the tarballs to the storage area of the cluster."""
    store = STORE[options.cluster].format(user=whoami())
    for tgz in tarball_paths():
        if not os.path.exists(tgz):
            if not DRY_RUN:
                sys.exit('%s not there yet, run with --tar first' % tgz)
            print('%s not there yet (--dry-run, --tar only printed above)' % tgz)
        target = store + '/' + os.path.basename(tgz)
        if options.cluster in XRDFS_RM:
            # this site cannot overwrite, so drop the old file (may not be there)
            run(['xrdfs', XRDFS_RM[options.cluster], 'rm', remote_path(target)],
                ignore_errors=True)
        run(['xrdcp', '-f', tgz, target])
        print('copied to ' + target)


def parse_args():
    parser = OptionParser()
    parser.add_option('-p', '--processor', help='processor name (also the hists/ sub directory)', dest='processor', default='')
    parser.add_option('-m', '--metadata', help='metadata name in metadata/<name>.json.gz', dest='metadata', default='')
    parser.add_option('-d', '--dataset', help='only these datasets, comma separated substrings', dest='dataset', default='')
    parser.add_option('-e', '--exclude', help='skip these datasets, comma separated substrings', dest='exclude', default='')
    parser.add_option('-f', '--file', help='take the job list from this file instead of scanning hists/', dest='file', default='')
    parser.add_option('-n', '--dry-run', help='write the list and the submit file but do not submit', action='store_true', dest='dry_run')
    parser.add_option('--site-config', help='JSON file with default values for the site options', dest='site_config', default='')
    # ---- site settings: the defaults come from SITE_DEFAULTS of the cluster ----
    parser.add_option('-c', '--cluster', help='cluster: %s (default $DECAF_CLUSTER or lpc)' % ', '.join(sorted(JDL)), dest='cluster', default=None)
    parser.add_option('--executable', help='job script to run (default run.sh)', dest='executable', default=None)
    parser.add_option('--cpus', help='request_cpus', dest='cpus', default=None)
    parser.add_option('--memory', help='request_memory in MB', dest='memory', default=None)
    parser.add_option('--disk', help='request_disk, e.g. 4000000 (default: let condor decide)', dest='disk', default=None)
    parser.add_option('--image', help='container image for +SingularityImage (lpc blob)', dest='image', default=None)
    parser.add_option('--accounting-group', help='accounting_group (kisti blob)', dest='accounting_group', default=None)
    parser.add_option('--proxy', help='x509 proxy file, "auto" to detect, "none" to skip (default auto)', dest='proxy', default=None)
    parser.add_option('--input-files', help='extra files to transfer with the job, comma separated', dest='input_files', default=None)
    parser.add_option('--exec-args', help='extra arguments appended to the executable, comma separated', dest='exec_args', default=None)
    parser.add_option('--extra', help='raw JDL line to append, repeatable', action='append', dest='extra', default=None)
    parser.add_option('--retries', help='max_retries for failing jobs (default 0, i.e. not set)', dest='retries', default=None)
    parser.add_option('--chunk', help='submit at most this many jobs per condor_submit (default 0 = all in one)', dest='chunk', default=None)
    parser.add_option('--submit-opts', help='extra options for condor_submit, e.g. "-spool"', dest='submit_opts', default=None)
    parser.add_option('--clean-logs', help='delete old condor logs of the resubmitted jobs', action='store_true', dest='clean_logs')
    # ---- staging of the environment the jobs unpack on the worker node ----
    parser.add_option('-t', '--tar', help='tar up decaf and the python site-packages', action='store_true', dest='tar')
    parser.add_option('-x', '--copy', help='xrdcp the tarballs to the storage of the cluster', action='store_true', dest='copy')
    options, args = parser.parse_args()

    if not options.processor:
        parser.error('-p/--processor is required')
    if not options.metadata and not options.file:
        parser.error('-m/--metadata is required (or give a job list with -f)')

    # command line wins over --site-config, --site-config wins over the defaults
    cfg = {}
    if options.site_config:
        with open(options.site_config) as fin:
            cfg = json.load(fin)
        unknown = [k for k in cfg if k not in SITE_KEYS]
        if unknown:
            sys.exit('unknown key(s) in %s: %s' % (options.site_config, ', '.join(sorted(unknown))))

    if options.cluster is None:
        options.cluster = cfg.get('cluster', os.environ.get('DECAF_CLUSTER', 'lpc'))
    if options.cluster not in JDL:
        sys.exit('unknown cluster "%s", pick one of: %s' % (options.cluster, ', '.join(sorted(JDL))))

    defaults = {
        'executable': 'run.sh',
        'disk': '',
        'proxy': 'auto',
        'input_files': '',
        'exec_args': '',
        'extra': [],
        'retries': '0',
        'chunk': '0',
        'submit_opts': '',
    }
    defaults.update(SITE_DEFAULTS[options.cluster])
    for key, value in cfg.items():
        if key != 'cluster':
            defaults[key] = value if key == 'extra' else str(value)
    for key, value in defaults.items():
        if getattr(options, key, None) in (None, ''):
            setattr(options, key, value)
    return options


def find_proxy(setting):
    """Return the proxy path to hand to condor, or '' if there is none."""
    if setting in ('none', 'no', ''):
        return ''
    if setting != 'auto':
        if not os.path.exists(setting):
            sys.exit('proxy file not found: ' + setting)
        return os.path.abspath(setting)
    candidates = [os.environ.get('X509_USER_PROXY', '')]
    try:
        candidates.append('/tmp/x509up_u%d' % os.getuid())
    except AttributeError:          # non posix
        pass
    for path in candidates:
        if path and os.path.exists(path):
            return os.path.abspath(path)
    print('WARNING: no x509 proxy found, submitting without one '
          '(use --proxy <file> or voms-proxy-init first)')
    return ''


def missing_jobs(options, outdir):
    """Jobs that have no output file in outdir yet."""
    if options.file:
        with open(options.file) as fin:
            jobs = [line.strip() for line in fin if line.strip() and not line.startswith('#')]
    else:
        with gzip.open('metadata/' + options.metadata + '.json.gz') as fin:
            jobs = list(json.load(fin).keys())

    # a job is done when outdir holds a file whose name without the extension
    # is the job name, whatever the extension is (.futures, .reduced, ...)
    done = set()
    if os.path.isdir(outdir):
        for name in os.listdir(outdir):
            done.add(name)
            done.add(os.path.splitext(name)[0])

    keep, wanted, skipped = [], options.dataset.split(','), options.exclude.split(',')
    for job in jobs:
        if options.dataset and not any(w in job for w in wanted if w):
            continue
        if options.exclude and any(s in job for s in skipped if s):
            continue
        if job in done:
            continue
        keep.append(job)
    return jobs, keep


def write_submit(options, path, listfile, outdir, proxy):
    inputs = [options.executable] + [f for f in options.input_files.split(',') if f]
    arguments = [options.metadata, '$(SAMPLE)', options.processor, options.cluster,
                 os.environ.get('USER', os.environ.get('LOGNAME', 'unknown'))]
    arguments += [a for a in options.exec_args.split(',') if a]

    proxy_lines = ''
    if proxy:
        proxy_lines = ('use_x509userproxy       = True\n'
                       'x509userproxy           = %s\n' % proxy)

    extra_lines = ''
    if str(options.disk):
        extra_lines += 'request_disk            = %s\n' % options.disk
    if int(options.retries) > 0:
        extra_lines += 'max_retries             = %s\n' % options.retries
    for line in (options.extra or []):
        extra_lines += line.rstrip() + '\n'

    jdl = JDL[options.cluster].format(
        sample='$(SAMPLE)',
        executable=options.executable,
        arguments=' '.join(arguments),
        input_files=', '.join(inputs),
        processor=options.processor,
        outdir=outdir,
        pwd=os.getcwd(),
        batchname='resubmit_' + options.processor,
        accounting_group=options.accounting_group,
        cpus=options.cpus,
        memory=options.memory,
        image=options.image,
        proxy_lines=proxy_lines,
        extra_lines=extra_lines,
        listfile=listfile,
    )
    with open(path, 'w') as fout:
        fout.write(jdl)


def clean_logs(processor, jobs):
    for sub in ('out', 'err', 'log'):
        folder = os.path.join('logs/condor/run', sub)
        if not os.path.isdir(folder):
            continue
        for name in os.listdir(folder):
            if not name.startswith(processor + '_'):
                continue
            if any(name.startswith('%s_%s_' % (processor, job)) for job in jobs):
                os.remove(os.path.join(folder, name))


def main():
    global DRY_RUN
    options = parse_args()
    DRY_RUN = bool(options.dry_run)
    if options.tar:
        make_tarballs()
    if options.copy:
        copy_tarballs(options)

    outdir = os.path.abspath(os.path.join('hists', options.processor))
    for folder in (outdir, 'logs/condor/run/out', 'logs/condor/run/err', 'logs/condor/run/log'):
        os.makedirs(folder, exist_ok=True)

    total, jobs = missing_jobs(options, outdir)
    tag = options.processor + ('_' + options.dataset.replace(',', '_') if options.dataset else '')
    listfile = 'missJobs_%s.txt' % tag
    with open(listfile, 'w') as fout:
        fout.write('\n'.join(jobs) + ('\n' if jobs else ''))

    print('%d jobs in %s' % (len(total), options.file or 'metadata/%s.json.gz' % options.metadata))
    print('%d outputs missing in %s' % (len(jobs), outdir))
    print('job list written to ' + listfile)
    if not jobs:
        print('nothing to resubmit')
        return 0

    proxy = find_proxy(options.proxy)
    submitfile = 'resubmit_%s.submit' % tag
    write_submit(options, submitfile, listfile, outdir, proxy)
    print('submit description written to %s (cluster: %s)' % (submitfile, options.cluster))

    if options.clean_logs:
        clean_logs(options.processor, jobs)

    if options.dry_run:
        print('--dry-run: nothing submitted, run "condor_submit %s" when happy' % submitfile)
        return 0
    if shutil.which('condor_submit') is None:
        print('condor_submit not found in PATH, nothing submitted')
        return 1

    chunk = int(options.chunk)
    chunks = [jobs] if chunk <= 0 else [jobs[i:i + chunk] for i in range(0, len(jobs), chunk)]
    rc = 0
    for idx, part in enumerate(chunks):
        if len(chunks) > 1:                       # one list file per condor_submit
            part_list = '%s.part%d' % (listfile, idx)
            with open(part_list, 'w') as fout:
                fout.write('\n'.join(part) + '\n')
            write_submit(options, submitfile, part_list, outdir, proxy)
            print('submitting %d/%d (%d jobs)' % (idx + 1, len(chunks), len(part)))
        cmd = ['condor_submit'] + [o for o in options.submit_opts.split() if o] + [submitfile]
        rc |= subprocess.call(cmd)
    return rc


if __name__ == '__main__':
    sys.exit(main())
