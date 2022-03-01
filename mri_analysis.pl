#!/usr/bin/perl
#This program is created by Rajan Paudel. The goal of this program is to compute nucleotide change in the MRI regions.
#This program can be excuted as perl mri_analysis.pl arg1
#arg1 can be  either of ag_r_85_; ct_r_85_;z_r_80;gc_r_77_;at_r_87_;gt_r_80_;ac_r_80_;
#For examples, ag_r_85 is typed to calculate the number of AG-rich regions with threshold more than 85%
#This program uses input files( the output files of regions.pl program) and also the VCF files(for each chromosome) from 1000 Genomes project
#This program outputs two files. The former file is used for further processing to give the final table.

$d = <@ARGV>;
$subs = substr $d,0,2;
$part1 = uc $subs;
$num = substr $d,5,2;
@arra = ("agct","ctag","gcat","atgc","gtac","acgt");
foreach $co (@arra){
	if ($co =~ /^$subs/){
		$ram = $co;
		}
	} 
$rai = uc $ram;
$part2 = substr $rai, 2,2;
@nu = split("", $rai);
$chromo = 1;$pairs ='';
for $ko(0..4){
	%{h.$ko} = 0;
	}
for $chromo(1..22){
	%ekpos = 0;
	$name = $d.$chromo;
	open(EK,"/home/rajan/5_1regions/$name")||die 'cannot open';
	while(<EK>){
		$_ =~ /chr\d+\t\d+\t(\d+)-(\d+)/;
		for $con($1..$2){
			$ekpos{$con}=1;
			}
	}
	close EK;
	open(IN,"zcat /home/patrick/2500GENOMES/ALL.chr$chromo.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz |")||die 'cannot find the file';
        while(<IN>){
                 if  ($_ =~ /^\d+\t(\d+)\t\w+\d+\t(\w)\t(\w).+;AF=(\d\.?\d*).+;AA=([ATGC])\|.+;VT=SNP/){
			$anc = $5;
			if ($anc eq $2){$der = $3;$dfr = $4;}
			elsif($anc eq $3){$der = $2;$dfr = 1-$4;}
			else{$der = 'X';$dfr = 2;}
			$pairs = $anc."\t".$der;
			if ($ekpos{$1}){
				$hall{$pairs}++;
        	                if($dfr > 0 && $dfr <= 0.20){$h0{$pairs}++;}else{;}
                	        if($dfr > 0.20 && $dfr <= 0.40){$h1{$pairs}++;}else{;}
                        	if($dfr > 0.40 && $dfr <= 0.60){$h2{$pairs}++;}else{;}
                       		if($dfr > 0.60 && $dfr <= 0.80){$h3{$pairs}++;}else{;}
                       		if($dfr > 0.80 && $dfr <= 1){$h4{$pairs}++;}else{;}
			}else{;}
		}
		print "$.\n";
	}
$. = 0;	
}	
close IN;
$namo = $d."table";
open(OUT,">$namo")|| die 'cannot create all';
foreach $am(keys %hall){
	print OUT "$am\t$h0{$am}\t$h1{$am}\t$h2{$am}\t$h3{$am}\t$h4{$am}\n";
	}
close OUT;

        open(IN1,"$namo")|| die 'cannot open file';
        while(<IN1>){
                chomp $_;
                @pi = split("\t",$_);
                if ($pi[0] =~ /$nu[0]/){
                                if($pi[1] =~ /$nu[1]/){
                                        for my $hi(2..6){
                                                ${mi.$hi} += $pi[$hi];
                                                }
                                        }
                                if($pi[1] =~ /$nu[2]/){
                                        for my $hi(2..6){
                                                ${mj.$hi} += $pi[$hi];
                                                }
                                        }
                                if($pi[1] =~ /$nu[3]/){
                                        for my $hi(2..6){
                                                ${mk.$hi} += $pi[$hi];
                                                }
                                        }
                                }
                if ($pi[0] =~ /$nu[1]/){
                                if ($pi[1] =~ /$nu[0]/){
                                        for my $hi(2..6){
                                                ${ni.$hi} += $pi[$hi];
                                                }
                                        }
                                if ($pi[1] =~ /$nu[2]/){
                                        for my $hi(2..6){
                                        ${nj.$hi} += $pi[$hi];
                                                }
                                        }
                                if ($pi[1] =~ /$nu[3]/){
                                        for my $hi(2..6){
                                        ${nk.$hi} += $pi[$hi];
                                                                               
                                                }
                                        }
                                }
                if ($pi[0] =~ /[$part1]/ and $pi[1] =~ /[$part2]/){
                        for my $hi(2..6){
                                ${ui.$hi} += $pi[$hi];
				}
			}
			
		if ($pi[0] =~ /[$part2]/ and $pi[1] =~ /[$part1]/){
                        for my $hi(2..6){
                                ${vi.$hi} += $pi[$hi];
				  }
                        }
                }
	close IN;
        for my $hi(2..6){
                push(@line1,${mi.$hi});
                push(@line2,${mj.$hi});
                push(@line3,${mk.$hi});
                push(@line4,${ni.$hi});
                push(@line5,${nj.$hi});
                push(@line6,${nk.$hi});
                push(@line7,${ui.$hi});
		push(@line8,${vi.$hi});
		${ti.$hi} = ${ui.$hi} / ${vi.$hi};
		${ti.$hi} = sprintf("%.2f",${ti.$hi});
		push(@line9,${ti.$hi});

        }
               open(OUT1,">Table_$part1");
		$" = "\t";
                $ranges = "0-20\t20-40\t40-60\t60-80\t80-100";
                print OUT1 "mutation\t$ranges\n$nu[0]>$nu[1]\t@line1\n$nu[0]>$nu[2]\t@line2\n$nu[0]>$nu[3]\t@line3\n\n$nu[1]>$nu[0]\t@line4\n$nu[1]>$nu[2]\t@line5\n$nu[1]>$nu[3]\t@line6\n\n$part1>$part2\t@line7\n$part2>$part1\t@line8\nR-value\t@line9\n";
                
       		 close OUT1;


