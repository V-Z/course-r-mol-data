#!/bin/bash

# After every Git commit affecting web...
# Use './web_update.sh all' to update also presentation and scripts and data

echo "Syncing 'rcourse' directory"
echo
rsync -arv --exclude=*~ --progress rcourse striga:~/soubory.trapa.cz/htdocs/ || { echo "Error!" && exit 1; }
echo

# Copy also presentation and scripts and data
if [ "$1" == "all" ]; then
	# Copy presentation
	echo "Syncing presentation"
	echo
	rsync -av --progress presentation/r_mol_data_phylogen.pdf striga:~/soubory.trapa.cz/htdocs/rcourse/ || { echo "Error!" && exit 1; }
	echo
	# Data
	# csv, html, in, nexus, nwk, r, txt
	fi

echo "Done!"

exit
