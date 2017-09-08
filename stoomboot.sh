echo "Launching " $1 " jobs of configuration " $2

for (( c=1; c<=$1; c++ ))
do  
   echo "cd ${PWD}
./dPDF ${c} ${2}" > ${c}.run
   qsub -q short ${c}.run
done
