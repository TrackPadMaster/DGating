#!/bin/bash

echo "Welcome to the first DGating script!"
echo "Brown's thesis only looks at p/q, and we follow suit here"
echo "I will need a P.txt and a Q.txt to work"


#starttime='date +%s'
# First we read the text files into arrays
readarray Pall < ./P.txt
readarray Qall < ./Q.txt
# q is just for simple indexing
echo ${Pall[@]}
echo ${Qall[@]}
q=0
# Cycle through every value of Pall
for p in ${Pall[@]}; do
	# For some infuriating reason, the readarray also saves spaces after the numbers
	# For most functions, this is just fine
	# For bc, the world is ending
	# A quick division by 1 (since they're integers) will fix this
	tempP=$(($p/1))
	tempQ=${Qall[$q]}
	tempQ=$(($tempQ/1))
	q=$(($q+1))
	# "scale" here decides how many digits the ratio will keep
	ratio=$(bc <<< "scale=5;($tempP/$tempQ)")
	# Now we have the ratio that we'll need to save
	awk -v ratio=$ratio '{if($3=="pqratio") {print ratio" !	pqratio"} else {print $0}}' GatingInStart > ./GatingIn0

	# GFortranRunner uses 40 phi points
	./GFortranRunner.sh
	./DataCollector.sh > ${tempP}over${tempQ}pqratio.dat

	echo "I just finished ${tempP} over ${tempQ}"
done
#endtime='date +%s'
echo "All done!"
#runtime=$((endtime-starttime))
#echo "I took ${runtime} seconds to run"

