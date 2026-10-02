#!/usr/bin/env python3
"""
Make a metadata json that holds only the jobs listed in a missJobs file.

The output is written next to the input as <ORIGINALNAME>_<PROCESSNAME>_resubmit.json.gz,
where PROCESSNAME is taken from missJobs_<PROCESSNAME>.txt, so it can be passed to
run.py / run_condor.py with -m <ORIGINALNAME>_<PROCESSNAME>_resubmit.

Examples
--------
./make_resubmit_json.py -m 2022_private_v2 -j missJobs_hadmonotop2022_0802.txt
./make_resubmit_json.py -m metadata/2022_private_v2.json.gz -j hadmonotop2022_0802
"""
import gzip
import json
import os
import sys
from optparse import OptionParser


def metadata_path(name):
    """Accept a bare metadata name (as run.py takes it) or a path to the file."""
    if os.path.exists(name):
        return name
    for cand in ('metadata/%s.json.gz' % name, 'metadata/%s.json' % name):
        if os.path.exists(cand):
            return cand
    sys.exit('metadata not found: %s' % name)


def missjobs_path(name):
    """Accept the missJobs file itself or just the PROCESSNAME part."""
    if os.path.exists(name):
        return name
    cand = 'missJobs_%s.txt' % name
    if os.path.exists(cand):
        return cand
    sys.exit('missJobs file not found: %s' % name)


def split_ext(path):
    base = os.path.basename(path)
    for ext in ('.json.gz', '.json'):
        if base.endswith(ext):
            return base[:-len(ext)], ext
    return os.path.splitext(base)


def main():
    parser = OptionParser()
    parser.add_option('-m', '--metadata', help='input json (name in metadata/ or a path)', dest='metadata', default='')
    parser.add_option('-j', '--jobs', help='missJobs file (or its PROCESSNAME)', dest='jobs', default='')
    parser.add_option('-o', '--outdir', help='output directory (default: same as the input json)', dest='outdir', default='')
    options, args = parser.parse_args()
    if not options.metadata or not options.jobs:
        parser.error('-m/--metadata and -j/--jobs are required')

    inpath = metadata_path(options.metadata)
    jobpath = missjobs_path(options.jobs)

    origname, ext = split_ext(inpath)
    process = os.path.splitext(os.path.basename(jobpath))[0]
    if process.startswith('missJobs_'):
        process = process[len('missJobs_'):]

    opener = gzip.open if inpath.endswith('.gz') else open
    with opener(inpath, 'rt') as fin:
        samplefiles = json.load(fin)

    with open(jobpath) as fin:
        jobs = [line.strip() for line in fin if line.strip() and not line.startswith('#')]

    # keep the order of the original json
    wanted = set(jobs)
    resubmit = {k: v for k, v in samplefiles.items() if k in wanted}
    unknown = [j for j in jobs if j not in samplefiles]

    outdir = options.outdir or os.path.dirname(inpath)
    outpath = os.path.join(outdir, '%s_%s_resubmit%s' % (origname, process, ext))
    with opener(outpath, 'wt') as fout:
        json.dump(resubmit, fout, indent=4)

    print('%d jobs in %s' % (len(jobs), jobpath))
    print('%d of %d entries of %s kept' % (len(resubmit), len(samplefiles), inpath))
    if unknown:
        print('WARNING: %d jobs not found in the input json:' % len(unknown))
        for j in unknown:
            print('  ' + j)
    print('written to %s' % outpath)


if __name__ == '__main__':
    main()
