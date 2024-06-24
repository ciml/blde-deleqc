#$ -S /bin/bash
#$ -q linux.q
#----------------------------------------------------------
# Parameters:
# 1- seed
# 2- F
# 3- genL
# 4- submission id (only used to create a folder)
#----------------------------------------------------------

# static parameters
cr=0.9
tol=0.0001
eps=0.0001

mkdir $4
cp *.cpp $4
cp *.h $4
cd $4

p=1
neval=1000

for varL in {1..4}; do

	for popL in 10 20 30 50 100; do
	
		for genF in 50 100 200 500; do

			# blde
			echo "blde"
			g++ blde.cpp funcoes.cpp -o cec2016
			
			for fobj in 3 4 9 25 26; do
				echo "./cec2016 -genL $3 -popL $popL -popF $popL -genF $genF -F $2 -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $1 | bzip2 -9 -c > cec2016_blde_${3}_${popL}_${genF}_${2}_${varL}_${fobj}_${1}.out.bz2"
				./cec2016 -genL $3 -popL $popL -popF $popL -genF $genF -F $2 -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $1 | bzip2 -9 -c > cec2016_blde_${3}_${popL}_${genF}_${2}_${varL}_${fobj}_${1}.out.bz2
			done
			
			# karla
			echo "karla"
			g++ karla.cpp funcoes.cpp -o cec2016
			
			for fobj in 20 21 22 23 24; do
				echo "./cec2016 -genL $3 -popL $popL -popF $popL -genF $genF -F $2 -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $1 | bzip2 -9 -c > cec2016_karla_${3}_${popL}_${genF}_${2}_${varL}_${fobj}_${1}.out.bz2"
				./cec2016 -genL $3 -popL $popL -popF $popL -genF $genF -F $2 -CR $cr -Var $varL -tol $tol -eps $eps -func $fobj -seed $1 | bzip2 -9 -c > cec2016_karla_${3}_${popL}_${genF}_${2}_${varL}_${fobj}_${1}.out.bz2
			done

		done

	done

done

cp *.out.bz2 ../results/

cd ..
