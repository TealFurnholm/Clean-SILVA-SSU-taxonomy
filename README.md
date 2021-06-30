# Clean-SILVA-SSU-taxonomy
Use my curated taxonomy database to get standardized lineages for the SILVA SSU genes
*make sure you've created the latest taxonomy database, described here: https://github.com/TealFurnholm/Universal-Taxonomy-Database

## perl FIX_SILVA.pl
This script: 
- downloads the current SILVA SSU database
- obtains the %identity and NCBI taxon identifier for each sequence
- removes duplicates and contained sequences @ 100% identity
- converts IUPAC >> "N" and "U" >> "T"
- creates the lowest common ancestor for each sequence using my curated taxonomy database
- sequence header given to the sequence with the highest scoring named organism, along with the %ID score, its Taxon ID, and its LCA
- sequence id mapping contains both the longest (absorbs lower LCAs) and shortest (worst case) LCA, and all of the taxon IDs

## inside the script
### auto-download SILVA's SSU gene metadata:
	- Download: https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/full_metadata/SILVA_138.1_SSURef.full_metadata.gz
	- includes the % identity of the genes annotation
	- includes the ncbi taxon identifier of the annotation
### auto-download SILVA gene sequences:
	- Download: https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz
### Get taxonomy, create LCAs: 
		- Convert any IUPAC code to N, convert U to T
			- comics reformat --fixjunk iupacton overwrite=t utot=t in=SILVA_138.1_SSURef_tax_silva.fasta out=SILVA_iupac.fasta;
		- Clean up the header
			- grep -oiP "^\S+" SILVA_iupac.fasta > SILVA_hedtrim.fasta;
		- Sort long >> short length
			- comics dedupe sort=length absorbrc=f absorbmatch=f absorbcontainment=f in=SILVA_hedtrim.fasta out=SILVA_sort.fasta;
		- Remove duplicates and contained sequences
			- comics dedupe sort=length absorbcontainment=t mergenames=t mergedelimiter=+ exact=f overwrite=t in=SILVA_sort.fasta out=SILVA_138.1_SSU.fasta;
			- Output of dedupe (version 138.1):
				Input: 		2224740 reads 		3183581141 bases.
				Duplicates: 	414296 reads (18.62%) 	603351762 bases (18.95%) 	0 collisions.
				Containments: 	157779 reads (7.09%) 	223893347 bases (7.03%) 	43685324995 collisions.
				Result: 	1652665 reads (74.29%) 	2356336032 bases (74.02%)
		- If multiple sequences/taxon identifiers - get the best (by % identity and by named species if available)
		- Create the shortest and longest LCA
		- Output sequences and mapping.
