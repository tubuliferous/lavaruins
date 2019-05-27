BASEDIR=$(dirname "$0")
cd $BASEDIR
nohup python lavaruins.py &
sleep 3s
open -a Firefox 'http://127.0.0.1:8050'