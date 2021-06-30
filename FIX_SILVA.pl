#use warnings;

### GET LATEST SILVA DATABASE ###

# SSU FASTA SEQS
qx{wget --spider --quiet -r --no-parent https://www.arb-silva.de/no_cache/download/archive/current/Exports/ -O tmp.txt};
$x = qx{cat tmp.txt};
$x =~ /fileadmin.silva_databases.current.Exports.(SILVA_([\d\.]+)_SSURef_tax_silva.fasta.gz)/i;
$file='https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/'.$1;
$version=$2;
qx{wget $file -O "silva_ssu.fasta.gz"};



### CLEAN SEQS ###
qx{comics reformat --fixjunk iupacton overwrite=t utot=t in=silva_ssu.fasta.gz out=silva_iupac.fasta};
$cmd = "grep -oiP \"^\\S+\" silva_iupac.fasta";
qx{$cmd > silva_hedtrim.fasta };
qx{comics dedupe sort=length absorbrc=f absorbmatch=f absorbcontainment=f overwrite=t in=silva_hedtrim.fasta out=silva_ssu_sort.fasta};
$outdb="SILVA_".$version."_dedup.fasta";
qx{comics dedupe sort=length absorbcontainment=t mergenames=t mergedelimiter=+ exact=f overwrite=t in=silva_ssu_sort.fasta out=$outdb};



# SSU METADATA
qx{wget --spider --quiet -r --no-parent https://www.arb-silva.de/no_cache/download/archive/current/Exports/full_metadata/ -O tmp.txt};
$x = qx{cat tmp.txt};
$x =~ /(SILVA_[\d\.]+_SSURef.full_metadata.gz)/i;
$file='https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/full_metadata/'.$1;
qx{wget $file -O "silva_ssu_meta.txt.gz"};
qx{gunzip silva_ssu_meta.txt.gz};



#use warnings;
$intx = "TAXONOMY_DB_2021.txt";
open(INTX,$intx)||die;
while(<INTX>){
	if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
	$tid=shift(@stuff);
	$lin=join(";",@stuff);
	$PHY{$tid}=$lin;
}


$inme = "silva_ssu_meta.txt";
open(INME,$inme)||die;
while(<INME>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_);
	$sid=$stuff[0];
        $tid=$stuff[$#stuff];
	$SID2TID{$sid}=$tid;
	$SID2SCO{$sid}=$stuff[43];
}


$/=">";
$outseq="SILVA_".$version."_SSU_NR100.fasta";
$outid="SILVA_".$version."_SSU_NR100_IDS.txt";
open(OUTSEQ,">",$outseq)||die "unable to open $outseq:$!\n";
open(OUTID,">",$outid)||die "unable to open $outid:$!\n";
print OUTID "Best_SilvaID\tBest_TaxonID\tShort_LCA\tLongest_LCA\tAll_TaxonIDs\n";
open(INFA, $outdb)||die "unable to open $outdb:$!\n";
while(<INFA>){
	if($_!~/\w/){next;}
	$_=uc($_);
	@stuff=split("\n",$_);
	$header=shift(@stuff);
	$seq=join("",@stuff);
	$seq=~s/[^A-Z]//g;

	@IDS=split('\+',$header);
	@PHYL=();	@ALT=();
	%TIDS=(); 
	$max=0;	  	$altmax=0;
	$best=''; 	$altbest='';
	$bid='';	$altbid='';
	foreach my $id (@IDS){
		$id=~/^([^\.]+)/;
		$sid=$1;
		$tid = $SID2TID{$sid};
		if(!exists($PHY{$tid})||$PHY{$tid}!~/\w/){ next;}
		$TIDS{$tid}=1;
		$lin = $PHY{$tid};
		$sco = $SID2SCO{$sid};
		push(@PHYL,$lin); 
		if($sco>$max){$max=$sco; $best=$sid; $bid=$tid;}
		if($lin=~/^(ARCHAEA|BACTERIA|MONA|EUKARYOTA)/){
			push(@ALT,$lin);
			if($sco>$altmax){$altmax=$sco; $altbest=$sid; $altbid=$tid;}
		}
	}

	#LOTS OF QUIDDAM (uncultured) IN SILVA
	#IF POSSIBLE USE NAMED ORGANISMS
	$ac=@ALT; if($ac>0){@PHYL=@ALT; $bid=$altbid; $max=$altmax; $best=$altbest;}

	#GET THE LONGEST and SHORTEST LCA
	@PTIPS=();
	@LINS=@PHYL;
	while($LINS[0]=~/[A-Z]/){
	        $lca = shift(@LINS);
	        if(grep(/$lca\;|$lca$/, @PTIPS) || grep(/$lca\;|$lca$/, @LINS)){ next; }
	        else{ push(@PTIPS,$lca);}
	}
	$shortlca = MakeLCA(@PHYL);
	$longlca  = MakeLCA(@PTIPS);
	print OUTSEQ ">".$best."|".$bid."|".$max."|".$shortlca."\n$seq\n";
	print OUTID "$best\t$bid\t$shortlca\t$longlca\t$bid";
	foreach my $tid (keys %TIDS){if($bid != $tid){print OUTID "|".$tid;}}
	print OUTID "\n";
}


sub MakeLCA{
        my @ARR=();
        %seen=();
        @array1 = @_;
        @array1 = grep { !$seen{$_}++ } @array1;

        #get the kingdoms, JIC lineage is NCA
        %LET=();
	$f=0;
        foreach my $lin (@array1){
		if($lin=~/HAZENIA/){$f=1;}
                if($lin !~ /^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/i){next;}
                $lin =~ /^(.)/; 
		$LET{$1}=1;
        }
        @LET=();
        foreach my $let (sort(keys %LET)){push(@LET,$let);}
        $let=join("",@LET);
	if($f==1){print "let $let array1 @array1\n";}

        $len1 = @array1;
        if($len1 == 1){$LCA = $array1[0]; }
        elsif($len1 > 1){
                $first = $array1[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @array1);
                        $len2 = @matched;
                        if($len2 == $len1){push(@ARR, $alevel);}
                        else{last;}
                }
                $len3 = @ARR;
                if($len3 > 1){$LCA = join(";", @ARR);}
                elsif($len3==1){$LCA = $ARR[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }
	if($f==1){print "let $let lca $LCA\n";}

        #add kingdoms to NCA
        if($LCA eq "NCA"){ $LCA.="-".$let; }
	if($f==1){print "let $let lca $LCA\n";}
        return($LCA);
}

