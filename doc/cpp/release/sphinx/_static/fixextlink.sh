#!/bin/bash

# get the list of lines that need to be fixed
IFS=$'\n'
matchlist=($(grep "href=\"http" _build/html/*.html))

# loop through those lines
for i in "${matchlist[@]}"; do

	# parse the line for the file name
	file=${i%%:*}
	
	# parse the line for the url
	urllist=($(echo $i | awk '$0=$2' FS=href\=\" RS=\>))

	# loop through urls
	for u in "${urllist[@]}"; do

		# remove last character
		umod=${u%?}
	
		# take care of / for sed
		umod=${umod//\//\\/}
		
		# is it a link to a NIST website 
		if [[ $u != *"nist.gov"* ]]; then
			
			# replace the text with the appropriate javascript command
			sed -ie "s/href=\\\"${umod}\\\"\>/href=javascript:redirectalert(\'${umod}\')\>/g" $file
		
		fi

	done

done
