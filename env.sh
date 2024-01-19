export PYTHONPATH=~/.local/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$PWD/analysis:$PYTHONPATH
export PYTHONWARNINGS="ignore"
#export PATH=~/.local/bin:$PATH
cat $HOME/private/$USER.txt | voms-proxy-init -voms cms --valid 140:00 -pwstdin
