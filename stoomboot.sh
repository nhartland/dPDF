echo "Launching " $1 " jobs of configuration " $2
xbase=${2##*/}
xpref=${xbase%.*}
for (( c=0; c<$1; c++ ))
do  
   echo "cd ${PWD}
./dPDF ${2} ${c}" > ${xpref}_${c}.run
   qsub -q short ${xpref}_${c}.run
done
