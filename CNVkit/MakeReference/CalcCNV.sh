#!/bin/bash

# Read from a file

# Obtains command line arguments
argsInput=$@


# Initializes options to default
outputPath="./"
referenceFile=""
fullGenomeComp=false
nextArgParam=false
rankingLength=0
shortMode=false


# If there is no user input, print out instructions
if [ $# -eq 0 ]
	then
	echo "Usage: CalcCNV.sh [options] [file1 [file2 ... ]]"
	echo "The following options are availible"
	echo " 	-l 	The number of files to display [default: Ranks all that are input]"
	echo " 	-r 	Reference file for bedtools intercet (Must be sorted)"
	exit
fi

# If the given path could not be found, quit the program
function checkFile() {

	# Checks file
	if [ $2 == "File" ]
		then

		# If the file cannot be read
		if [ ! -r $1 ]
			then

			# If the file does not exist, print an error and quit
			if [ ! -e $1 ]
				then
				echo "ERROR: The file $1 could not be found. Skipping..."
				return 1
				
			# If the file exists but is not readable, print an error and quit
			else
				echo "ERROR: No read permissions for $1. Skipping..."
				return 1
				
			fi
		fi

		# Checks Directory
	elif [ $2 == "Dir" ]
		then

		if [ ! -d $1 ]
			then
			echo "ERROR: The path $1 Could not be found"
			return 1
			
		fi
	fi

	return 0
}

# Calculates how significantly the copy number deviates from standard diploid
function calcDeviation() {

	rawPloidy=$1
	echo $1
	totalOffset=0
	# Obtains Ploidy from file
	for linePloidy in ${rawPloidy[*]}
	do

		let varOffNormal=$linePloidy-2
		# Obtains absolute value
		if [ $varOffNormal -lt 0 ]
			then
			let varOffNormal=$varOffNormal*-1
		fi
		let totalOffset=$totalOffset+$varOffNormal

	done
	echo $totalOffset

}

# Creates a string with the desired number of elements
function initializeArray() {

	length=$1
	filler=$2
	for ((i=0;i<$length;i++))
	do
		arrayString="$arrayString $2"
	done
	echo $arrayString

}

function numOfFiles() { 

	echo $# 

}


function moveUp() {

	for ((i=${#lowestPanelVar[*]}-2;i>=$1;i--))
	do
		echo "Moving up"
		let nextIndices=$i+1
		lowestPanelVar[$nextIndices]=${lowestPanelVar[$i]}
		lowestFiles[$nextIndices]=${lowestFiles[$i]}

	done

}


# Treats each argument in the file
currentArg=-1
nextArgTreated=false
argsInput=($argsInput)
for arg in ${argsInput[*]}
do

	# The parameter will be treated as an input file name
	((currentArg++))
	if [ $nextArgTreated == true ]
	then
		nextArgTreated=false

	elif [ $arg == "-r" ] && [ $((currentArg+1)) -le ${#argsInput[*]} ]
	then
		checkFile ${argsInput[$((currentArg+1))]} "File"
		referenceFile=${argsInput[$((currentArg+1))]}
		nextArgTreated=true

	elif [ $arg == "-l" ] && [ $((currentArg+1)) -le ${#argsInput[*]} ]
	then
		rankingLength=${argsInput[$((currentArg+1))]}
		nextArgTreated=true
	elif [ $arg == "-s" ]
	then
		shortMode=true
	else
		checkFile $arg "File"
		if [ $? -eq 0 ] 
		then
			inputFiles="$inputFiles $arg"
		fi

	fi

done

# Quits if no reference file was provided
if [[ $referenceFile == "" ]]
	then
	echo "ERROR: No reference File was provided"
	echo "This is a file containing annotated regions"
	echo "Specify a referenece file using '-r' "
	exit
fi

# Converts strint to array
inputFiles=($inputFiles)
if [ $rankingLength -eq 0 ] 
then
	rankingLength=${#inputFiles[*]}
fi

lowestPanelVar=( $(initializeArray $rankingLength 100000) )
lowestTotalVar=( $(initializeArray $rankingLength 100000) )
lowestFiles=( $(initializeArray $rankingLength None) )

# The working bit of the program
for file in ${inputFiles[*]}
do

	# Makes a bed file from the input file, with the 4th column as the ploidy of each region
	# Removes the header, and uses bedtools intersect to only obtain the panneled regions
	# The 4th column (containing the ploidy) is then cut from the input, and stored as a list
	filePanelPloidy=$(cut -f 1,2,3,4 $file | tail -n +2 | bedtools intersect -b $referenceFile -a stdin | cut -f 4)
	
	panelOffset=0
	# Calculates ploidy deviation
	for ploidy in ${filePanelPloidy[*]}
	do

		# Normalizes ploidy against dipliod organisms
		let varOffNormal=$ploidy-2
		# Obtains absolute value

		if [ $varOffNormal -lt 0 ]
			then
			let varOffNormal=$varOffNormal*-1
		fi
		let panelOffset=$panelOffset+$varOffNormal

	done

	totalOffset=0
	fileTotalPloidy=$(cut -f 4 $file | tail -n +2)
	for ploidy in ${fileTotalPloidy[*]}
	do

		# Normalizes ploidy against dipliod organisms
		let varOffNormal=$ploidy-2
		# Obtains absolute value

		if [ $varOffNormal -lt 0 ]
			then
			let varOffNormal=$varOffNormal*-1
		fi
		let totalOffset=$totalOffset+$varOffNormal

	done

	# Compares the current offset to the highest-ranking versions stored
	# Shifts the list as needed to order the new file apropriately
	for ((i=0;i<${#lowestPanelVar[@]};i++))
	do
		# Compares CNVs in panel regions. If the offset is the pannel region is less than the current entry, it is added to the list
		if [ $panelOffset -le ${lowestPanelVar[$i]} ]
		then

			
			while [ $i -lt $((${#lowestPanelVar[@]}-1)) ] && [ $panelOffset -eq ${lowestPanelVar[$i]} ] && [ $totalOffset -gt ${lowestTotalVar[$i]} ]
			do
				((i++))
				
			done

			# Shifts the array upwards as nessisary
			for ((j=${#lowestPanelVar[*]}-2;j>=$i;j--))
			do
				let nextIndices=$j+1
				lowestPanelVar[$nextIndices]=${lowestPanelVar[$j]}
				lowestTotalVar[$nextIndices]=${lowestTotalVar[$j]}
				lowestFiles[$nextIndices]=${lowestFiles[$j]}

			done
			#Sets the new variables
			lowestPanelVar[$i]=$panelOffset
			lowestTotalVar[$i]=$totalOffset
			lowestFiles[$i]=$file
			break

		fi
			
	done

done

# Prints out the ranked order of each input file
if [ $shortMode == true ]
	then 
	for ((i=0;i<${#lowestFiles[@]};i++))
	do
		noExtName=${lowestFiles[$i]##*/}
		baseName=${noExtName%%.*}
		echo -n "$baseName "
	done
	echo ""
else
	echo ""
	echo "The files, ranked in order are"

	for ((i=0;i<${#lowestFiles[@]};i++))
	do
		echo "#$((i+1)) : ${lowestFiles[$i]} with a score of ${lowestPanelVar[$i]}:${lowestTotalVar[$i]}"
	done

fi
