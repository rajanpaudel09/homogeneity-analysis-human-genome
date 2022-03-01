#!/usr/bin/perl
#This program is created by Rajan Paudel. The goal of this program is to calculate different MRI regions inside the human genome. It also computes differnet characterstics of the MRI regions.
#This file uses sequences from UCSC genome browser as input files
# The program can be run as, perl regions.pl arg1,   where arg1 is can be 'AT_87,GC_77,AG_85,CT_85,AC_80,GT_80'
# The input file is the fasta file from hg19 reference genome from UCSC genome browser
# There are two types of output files. One is summary file containing information about the MRI region.
# The other type of files are the files containing MRI regions, their co-ordiantes and length for each chromosomes

$d = <@ARGV>;
@ele = split("",$d);
$str1 = substr $d,0,2;
$str2 = substr $d,3,2;
$into = int $str2;
$int = $into/100;

open(OUT,">Summary_$d");
$ch = 1;
for $ch(1..22){
	#system("tail -n +2 chr$ch.fa | paste -s -d \"\" >chr$ch");
	$name = $str1."_r_".$str2."_".$ch;
	open(OUT1,">$name")||die 'cannot do open';
	open(IN,"/home/rajan/the1120/chr/chr$ch")|| die "Couldn't open file: $!";
	$ini = '';$tot = 0;$p = 0;$n=1;$se = '';$seq ='';
	while (<IN>) {
		$whole = $_;
		$leng = () = $_ =~ /[GgCcAaTtN]/g;
	}
	print "$leng\n";
	close IN;
	while(){last if ($n > $leng);
		$whole =~ /\G(.{100})/gc;$ini = $1;$n = $n+100;
		$na = () = $ini =~ /[$str1]/gi;
		if ($na > $into){$m = $n-100;
               		$seq = $ini;
			START:
			$whole =~ /\G(.{10})/gc;
        		$str = $1;$pseq = $seq;$pn = $n;$nall = $pn - $m;
			$seq = $seq.$str;
        		$nam = () = $seq =~ /[$str1]/gi;$n = $n + 10;$nat = $nall +10;$at = $nam/$nat;
        			if ($at > $int){goto START;
				}else{$p = $p+$pn-$m;$wseq = $wseq.$pseq; print OUT1 "chr$ch\t$nall\t$m-$pn\t$pseq\n";$count++;
			}
		}
	}
	$mo = $p/$n;
	close OUT1;
print "$.\n";
}
my $all = () = $wseq =~ /[agtcAGCTN]/g;
my $no1 = () = $wseq =~ /[$ele[0]]/ig;
my $no2 = () = $wseq =~ /[$ele[1]]/ig;
my $rep = () = $wseq =~ /[agtc]/g;
my $avg = $all/$count;
my $rat = $rep/$all*100;
my $per1 = $no1/$all*100;
my $per2 = $no2/$all*100;
print "all -$all\tno1 $no1\t no2 $no2\t rep $rep\t avg $avg\t rat $rat\t %1 $per1\t%2 $per2\n";
$avg = sprintf("%.0f",$avg);
$rat = sprintf("%.2f",$rat);
$per1 = sprintf("%.2f",$per1);
$per2 = sprintf("%.2f",$per2);
print OUT "Total no of regions = $count\n Total lenth of $str1 -rich MRI region = $all\n% of $ele[0] = $per1 \n % of $ele[1] = $per2 \n Average length of $d  = $avg\n";
close OUT;



