#!/usr/bin/perl
#Program created by Rajan Paudel. Its goal is to calculate the number of SNPs with ancestral and mutant allele information in the human genome using data from 1000 Genomes project
#Command line for starting program: perl wg_analysis.pl
#The input files are the VCF files for each chromosome from 1000 genomes project
#There are two output files 'whole_genome' and 'whole_genome_NR_values'.
#The output files 'whole_genome' which gives summary of ancestrl and derived allele change
#The output files 'whole_genome_NR_values' gives the summary of 'N' and 'R' values for each bin

for my $chromo(1..22){
        open(OUT,">whole_genome");
        open(IN,"zcat /home/patrick/2500GENOMES/ALL.chr$chromo.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz |")||die 'cannot find the file';
        while(<IN>){
                if  ($_ =~ /^\d+\t(\d+)\t\w+\d+\t(\w)\t(\w).+;AF=(\d\.?\d*).+;AA=([ATGC])\|.+;VT=SNP/){
			$anc = $5;
        		if ($anc eq $2){$der = $3; $dfr = $4;}
        		elsif($anc eq $3){$der = $2;$dfr = 1-$4;}
        		else{$der = 'X';$dfr = 2;}
        		my $pairs = $anc."\t".$der;
        		$hall{$pairs}++;
        		if($dfr > 0 && $dfr <= 0.20){$h1{$pairs}++;}else{;}
        		if($dfr > 0.20 && $dfr <= 0.40){$h2{$pairs}++;}else{;}
        		if($dfr > 0.40 && $dfr <= 0.60){$h3{$pairs}++;}else{;}
        		if($dfr > 0.60 && $dfr <= 0.80){$h4{$pairs}++;}else{;}
        		if($dfr > 0.80 && $dfr <= 1){$h5{$pairs}++;}else{;}
			}
		else{;}
		}
	}
foreach $am(keys %hall){
        print OUT "$am\t$h1{$am}\t$h2{$am}\t$h3{$am}\t$h4{$am}\t$h5{$am}\n";
        }
close IN;
close OUT;
open(OUT1,">whole_genome_NR_values");
open(DAT,"whole_genome")|| die 'cannot open file';
        while(<DAT>){
                chomp $_;
		@a = split ("\t", $_);
		if($a[0] =~ /[GC]/ && $a[1] =~ /[AT]/){
                       for $hi(2..6){
                       		${mi.$hi} += $a[$hi];
			}
		}
		if($a[0] =~ /[AT]/ &&   $a[1] =~ /[GC]/){
                       for $hi(2..6){
                       		${ni.$hi} += $a[$hi];
			}
		}
	}
close DAT;
for $hi (2..6){
	push(@line1,${mi.$hi});
	push(@line2,${ni.$hi});
	${r.$hi} = ${mi.$hi}/${ni.$hi};
	${n.$hi} = ${r.$hi}*58/42;
	${r.$hi} = sprintf("%.2f",${r.$hi});
	${n.$hi} = sprintf("%.2f",${n.$hi});
	push(@line3,${r.$hi});
	push(@line4,${n.$hi});
	}
$" = "\t";
$ranges = "0-20\t20-40\t40-60\t60-80\t80-100";
print OUT1 "\t$ranges\nGC>AT\t@line1\nAT>GC\t@line2\nR-value\t@line3\nN-value\t@line4\n";
close OUT1;







