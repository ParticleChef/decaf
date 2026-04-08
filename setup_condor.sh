source /eos/user/t/taiwoo/miniconda3/etc/profile.d/conda.sh
conda activate py38
voms-proxy-init -voms cms -valid 24:00
voms-proxy-init -voms cms -valid 24:00   -out /eos/user/t/taiwoo/decaf/analysis/proxy/x509up_u$(id -u)

tar -czf analysis.tgz \
    --exclude='*.root' \
    --exclude='*.log' \
    --exclude='*.out' \
    --exclude='*.txt' \
    --exclude='*.err' \
    --exclude='*.pdf' \
    --exclude='*.png' \
    --exclude='*.pyc' \
    --exclude='*.futures' \
    --exclude='__pycache__' \
    --exclude='condor_out' \
    --exclude='*.reduced' \
    --exclude='plots*' \
    analysis/

mv analysis.tgz condor/

module load lxbatch/eossubmit