dirs=$(ls -d PHI_*)

echo "# phi ; <v> ; error <v> ; diff ; error diff ; Peclet"

for dir in $dirs ; do

cat $dir/output.dat | awk '{if($1=="phi") {phi=$3} if($1=="<v>") {v=$3 ; getline ; dv=$2 ; getline ; d=$3 ; getline ; dd=$2 ; getline ; pe=$3 }}END{print phi"   "v"   "dv"    "d"     "dd"    "pe}' 


done
