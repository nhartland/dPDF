echo "Launching " $1 " jobs of configuration " $2

for (( c=0; c<$1; c++ ))
do  
   echo "cd ${PWD}
./dPDF ${2} ${c}" > ${c}.run
   qsub -q short ${c}.run
done
