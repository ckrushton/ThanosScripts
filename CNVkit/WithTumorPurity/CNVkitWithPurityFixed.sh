#!/bin/bash

# This pipeline runs multiple samples through CNVkit if a different purity is presented for each file
# Requires a reference to already be constructed. See CNVkit instructions on how to build a reference


# Lists out the arguments required to run CNVkit
samplePurityFile="None"
reference="None"
scatter=false
heatmap=false
diagram=false
inputBAMs=()
inputVCFs=()
outputPath="None"

function runCall() {

	cnsFiles=$1
	cnrFiles=$2

	for bamFile in ${inputBAMs[*]}
	do

		cnvkitCallCom=""
	
		# Obtains the base name of the input
		bamNoExtensions=${bamFile%%.*}
		baseBAMName=${bamNoExtensions##*/}

		bamCns="${outputPath}${baseBAMName}.cns"
		bamCnr="${outputPath}${baseBAMName}.cnr"
		bamCallOutput="${outputPath}${baseBAMName}.call.cns"

		matchedVCF="None"
		#Pairs each .bam with a .vcf, if applicable
		if [ ${#inputVCFs[*]} -gt 0 ]
		then

			for vcfFile in ${inputVCFs[*]}
			do

				vcfNoExtension=${vcfFile%%.*}
				baseVCFName=${vcfNoExtension##*/}

				if [ $baseVCFName == $baseBAMName ]
				then

					matchedVCF=$vcfFile
					break
				
				fi

			done

			if [ matchedVCF != "None" ]
			then
				cnvkitCallCom+=" -v $matchedVCF"
			else
				echo "ERROR: No .vcf file was found for $bamFile. Proceeding..."

			fi

		fi

		samplePurity=-1
		bamNoUnderscore=${baseBAMName%_*}${baseBAMName#*_}
		#Finds a matching purity file
		while read purityLine;
		do
			purityTokens=($purityLine)
			purityTokenName=${purityTokens[0]}
			if [ $purityTokenName == $bamNoUnderscore ]
			then
				samplePurity=$(echo "1 - ${purityTokens[1]}" | bc -l )
				break
			fi
			
		done<$samplePurityFile

		# If sample purity was not found in the file, the .bam file is skipped
		if [ $samplePurity == -1 ]
		then
			echo "ERROR: No purity value could be found for $bamFile. Skipping..."
			continue
		fi

		cnvkitCallCom+=" -m clonal --purity $samplePurity"

		cnvkitCallCom+=" -o $bamCallOutput"
		cnvkit.py call $bamCns $cnvkitCallCom

		# Runs visualization commands on corrected .cns
		if [ $scatter == true ] 
		then 
			if [ $matchedVCF != "None" ]
			then

				cnvkit.py scatter $bamCnr -s $bamCallOutput -v $matchedVCF -o ${outputPath}${baseBAMName}.call.scatter.pdf
			else

				cnvkit.py scatter $bamCnr -s $bamCallOutput -o ${outputPath}${baseBAMName}.call.scatter.pdf
			fi
		fi

		if [ $heatmap == true ]
		then
		
			cnvkit.py heatmap $1 -o Sample-Wide-Heatmap.cns.pdf
			cnvkit.py heatmap $2 -o Sample-Wide-Heatmap.cnr.pdf
		fi

		if [ $diagram == true ]
		then

			cnvkit.py diagram $bamCnr -s $bamCallOutput -o ${outputPath}${baseBAMName}.call.diagram.pdf
		fi

	done
}



cmdArgs=($@)
# Processes command line arguments
currentArg=-1
argAlreadyTreated=false
for arg in ${cmdArgs[*]}
do

	((currentArg++))
	if [ $argAlreadyTreated == true ]
	then
		argAlreadyTreated=false
		continue
	fi

	# If an option is specified
	if [ ${arg:0:1} == "-" ]
	then

		# If an option is specified that requires an additional parameter
		if [ $arg == "-p" ] || [ $arg == "-r" ] || [ $arg == "-t" ] || [ $arg == "-a" ] || [ $arg == "-o" ]
		then

			# Checks to make sure there is a another argument
			if [ $((currentArg+1)) -lt ${#cmdArgs[*]} ]
			then

				# If the option specifies a path, skip the file exists check
				if [ $arg == "-o" ]
				then

					outputPath=${cmdArgs[$((currentArg+1))]}
					argAlreadyTreated=true
					continue
				fi

				# If the purity file is specified, set the variable
				if [ $arg == "-p" ]
				then 

					# Checks to make sure the file exists
					if [ ! -e ${cmdArgs[$((currentArg+1))]} ]
					then

						echo "The file $arg could not be found. Skipping"
					else					
						samplePurityFile=${cmdArgs[$((currentArg+1))]}	
					fi

				# If the reference file is specified, set it
				elif [ $arg == "-r" ]
				then

					# Checks to make sure the file exists
					if [ ! -e ${cmdArgs[$((currentArg+1))]} ]
					then

						echo "The file $arg could not be found. Skipping"
					else

						reference=${cmdArgs[$((currentArg+1))]}
					fi

				fi

				argAlreadyTreated=true
			# If there is no next argument for these parameters, return an error
			else

				echo "ERROR: No parameter was specified for $arg"
				exit
			fi

		# Sets the output options
		elif [ $arg == "--scatter" ]
		then

			scatter=true

		elif  [ $arg == "--heatmap" ]
		then

			heatmap=true

		elif [ $arg == "--diagram" ]
		then

			diagram=true

		# If an unsuported option was specified, print out an error message
		else

			echo "ERROR: The option $arg is not supported"
			exit
		fi

	# If there is no "-" at the start, the input is sorted into the apropriate array
	else
		
		# Checks to make sure the file exists
		if [ ! -e $arg ]
		then

			echo "The file $arg could not be found. Skipping"

		else

			fileExtension=${arg##*.}

			# Adds a .bam file to the list
			if [ $fileExtension == "bam" ]
			then
				inputBAMs+=($arg)

			elif [ $fileExtension == "vcf" ]
			then
				inputVCFs+=($arg)
			else 

				echo "$arg is an unknown file type, ignoring..."
			fi
		fi
	fi
done


# Checks the input files to ensure that all required files are specified
# If no .BAM files were provided, prints and error message and exits
if [ ${#inputBAMs[*]} -eq 0 ]
then
	echo "ERROR: No input .BAM(s) were provided"
	exit
fi

# Checks to ensure the reference file was specified and is valid
if [ $reference == "None" ]
then

	echo "ERROR: No reference file was provided"
	exit

elif [ ${reference##*.} != "cnn" ]
then

	echo "ERROR: The specified reference file may be invalid."
	echo "The file extension is ${reference##*.}, but \".cnn\" was provided"
	exit

fi

# Ensures a purity file was provided
if [ $samplePurityFile == "None" ]
then 
	echo "ERROR: No purity file was provided"
	exit
fi


# If there is default output path, set it 
if [ $outputPath == "None" ]
then

	chosenBAM="${inputBAMs[0]}"

	# If the input files are in the local directory, set the output to the local directory
	if [ $chosenBAM == ${chosenBAM%/*} ]
	then
		outputPath="./"

	else
		outputPath="${chosenBAM%/*}/"
	fi

fi

# Starts creating the command string for batch
cnvkitBatchCommands=" -r $reference --output-dir $outputPath -p 8"

source activate python2

# cnvkit.py batch ${inputBAMs[*]} $cnvkitBatchCommands 

# Creates a string with additional parmeters to pass to call

cnvCallCom=""

if [ $scatter == true ] 
then 
	cnvkitCallCom+=" --scatter"
fi

if [ $heatmap == true ]
then
	
	cnvkitCallCom+=" --heatmap"
fi

if [ $diagram == true ]
then

	cnvkitCallCom+=" --diagram"
fi


runCall ${outputPath}*.cns ${outputPath}*.cnr

source deactivate
exit
