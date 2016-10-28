#!/bin/bash
# Runs cnvkit.py call on a range of purities

source activate python2
cmdArgs=$@
for cnsFile in $cmdArgs
do

	cnsNoExtensions=${cnsFile%%.*}
	baseCnsName=${cnsNoExtensions##*/}
	mkdir baseCnsName

	for ((j=10;j<100;j=j+5))
	do
		purity="0.${j}"
		cnvkit.py call $cnsFile -m clonal --purity $purity -o ${baseCnsName}/${baseCnsName}.callP${j}.cns
		cnvkit.py scatter /home/crushton/CNVkit/TumorRefWalkthrough/Ref_From_Tumor_47/${baseCnsName}.cnr -s ${baseCnsName}/${baseCnsName}.callP${j}.cns -o ${baseCnsName}/${baseCnsName}-scatter.pdf
		cnvkit.py diagram /home/crushton/CNVkit/TumorRefWalkthrough/Ref_From_Tumor_47/${baseCnsName}.cnr -s ${baseCnsName}/${baseCnsName}.callP${j}.cns -o ${baseCnsName}/${baseCnsName}-diagram.pdf

	done
done

source deactivate
exit