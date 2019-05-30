#!/bin/bash
count_pngs=`ls -1 simulations/*.png 2>/dev/null | wc -l`
if [ $count_pngs != 0 ]
then
rm ./simulations/*png
fi
count_txts=`ls -1 simulations/*.txt 2>/dev/null | wc -l`
if [ $count_txts != 0 ]
then
rm ./simulations/*txt
fi
python ./src/main.py $1 $2 $3 $4
declare pids
# for pid in $(pgrep 'draw'); do
#     wait "$pid"
# done
# declare -a $(pgrep 'draw')
# pwait 'draw' --name
# wait pgrep -f draw
# start=0
# end=126
# while [ $end -le $2 ]; do
#     ./draw $start $end
#     (( start+=126 ))
#     (( end+=126 ))
# done
# echo $end
# ./draw $start $2
# if [ "$(ls -A simulations/plane_*png)" ]; then
#if ls ${simulations}/*.png;
count_pngs=`ls -1 simulations/*.png 2>/dev/null | wc -l`
if [ $count_pngs != 0 ]
then
  yes | ffmpeg -i simulations/plane_%03d.png simulations/simul.gif
fi
