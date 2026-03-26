voms-proxy-init -voms cms -valid 24:00   -out /eos/user/t/taiwoo/decaf/analysis/proxy/x509up_u$(id -u)

tar -czf analysis.tgz \
    --exclude='*.root' \
    --exclude='*.log' \
    --exclude='*.pdf' \
    --exclude='*.png' \
    --exclude='*.pyc' \
    --exclude='__pycache__' \
    --exclude='condor_out' \
    --exclude='hists' \
    --exclude='plots' \
    analysis/

mv analysis.tgz condor/

module load lxbatch/eossubmit