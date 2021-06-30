# Clean-SILVA-SSU-taxonomy
Use my curated taxonomy database to get standardized lineages for the SILVA SSU genes

* make sure you've created the latest taxonomy database, described here: https://github.com/TealFurnholm/Universal-Taxonomy-Database
* Once you have the taxonomy database, download "FIX_SILVA.pl" and run "perl FIX_SILVA.pl" on the shell command line
* OUTPUTS: 
	1. SILVA_[version]_SSU_NR100.fasta
	2. SILVA_[version]_SSU_NR100_IDS.txt

## about FIX_SILVA.pl
This script: 
- downloads the current SILVA SSU database
- obtains the %identity and NCBI taxon identifier for each sequence
- removes duplicates and contained sequences @ 100% identity
- converts IUPAC >> "N" and "U" >> "T"
- creates the lowest common ancestor for each sequence using my curated taxonomy database
- sequence header given to the sequence with the highest scoring named organism, along with the %ID score, its Taxon ID, and its LCA
- sequence id mapping contains both the longest (absorbs lower LCAs) and shortest (worst case) LCA, and all of the taxon IDs
*Note:
	- Output of dedupe (version 138.1):
		Input: 		2224740 reads 		3183581141 bases.
		Duplicates: 	414296 reads (18.62%) 	603351762 bases (18.95%) 	0 collisions.
		Containments: 	157779 reads (7.09%) 	223893347 bases (7.03%) 	43685324995 collisions.
		Result: 	1652665 reads (74.29%) 	2356336032 bases (74.02%)
	- If multiple sequences/taxon identifiers - get the best (by % identity and by named species if available)
	- Since the provenance of each sequences identification is unknown, I stuck with the short LCA, so if there are 3 identifiers: species, genus, and kingdom, the LCA will be to the kingdom only. There is no good balance here, as the species may share the same genus and the stronger evidence may be to call the LCA at the genus level. The longest LCA absorbs any lower LCAs, and then finds the common LCA.
