#! /usr/bin/env perl
# Count-STM - Programed by Mauricio J Lozano
use strict;
use warnings;

#variables
my %cod;
my %aa;
my @codons_a;
my @seqcodons;
my $tempfile;
my $counts;
my $head;

#vars used on subs
my %res_counts;
my %res_freqs;
my %res_rscus;

# abre archivos de barcodes y pools
open (CODONS, "<codons2.txt");
@codons_a = <CODONS>;
close (CODONS);

#indexa los codones, el número de cada codones por aa etc
foreach (@codons_a) {
	(my $aas, my $ncod, my $codon) = split(",", $_);
	chomp($codon);
	$cod{$aas}=$codon;
	$aa{$aas}=$ncod;
}
$head = "Genes/Codons";
foreach (sort keys %cod){
	$head = $head." ".$_;
}


#archivos
my $file1;
if (@ARGV){
	$file1 = $ARGV[0];
    } 
else{
    print "Please type in the file name: ";
    $file1 = <STDIN>;
}

$file1 =~ /(^.+)\./;
my $name = $1;
print "File name: ".$name."\n";
my $file2 = $name.".counts";
my $file3 = $name.".freqs";
my $file4 = $name.".rscu";

#input.dat and temp.file
$tempfile = readfile ($file1,"seq.temp");

#calculate counts
$counts = counts ($tempfile,$file2, $file3, $file4);


unlink "seq.temp";



#**********************************************
sub readfile {
	my @files=@_;
	
	my $input=$files[0];
    my $temp=$files[1];
	my @tempo;
	my $i=0;

	open (TEMP,">$temp");
	open (DATOS, "<$input");
	chomp(my @lines = <DATOS>);

	 foreach my $line (@lines) {
		if ($line =~ /^>/){
			$line =~ s/\s/_/g;
			$line =~ s/\'/_/g;
			push @tempo, "\n".$line."\n";
		}
		elsif($line =~ /^$/){
		}
		else {
			$line =~ s/\r//g;
            $line =~ s/-*$//;
			$line =~ tr/agtc/AGTC/;
			push @tempo, $line;
		}
	}
	$tempo[0] =~ s/^\n//;
	foreach my $line (@tempo){
		print TEMP $line;
	}
	close TEMP;
	close DATOS;
	return $temp;
}

#**********************************************
sub counts {
	my @files=@_;
	
	my $input=$files[0];
    my $count=$files[1];
    my $freq=$files[2];
    my $rscu=$files[3];

	open (COUNT,">$count");
	open (FREQCU,">$freq");
	open (RSCU,">$rscu");
	open (DATOS, "<$input");

	$res_counts{0} = $head;
    $res_freqs{0} = $head;
    $res_rscus{0} = $head;
    chomp(my @lines = <DATOS>);
    my $i=0;
	my $id;

    #calculate counts
	foreach my $line (@lines){
		my %codoncount;
		my %codoncount2;
		my @seqcount;

		if ($line =~ /^>/){
			$line =~ s/>//;
			$id = $line;
	        $i++;
		}
		else {
			foreach my $acod (sort keys %cod){
				$codoncount{$acod} = 0;
				$codoncount2{$cod{$acod}}=0;
			}

			@seqcodons = $line =~ /(...)/g;		
			foreach my $codon (@seqcodons){
				$codoncount2{$codon}++;
			}
			
			foreach (sort keys %cod){
				my $acod=$_;
				$codoncount{$acod} = $codoncount2{$cod{$acod}};
			}

			foreach my $keys (sort keys %codoncount){
				push @seqcount,$codoncount{$keys};
			}

			my $counted = join(' ',@seqcount);
			$res_counts{$i} = $id." ".$counted;

			#calculate freqs and rscus
			
			my @seqfreq;
			my @seqrscu;

			my %codfreqs = %codoncount;
			my %ncodons;
			my %sums;
			foreach my $key (sort keys %codfreqs){
				$key =~ /.{3}(\w\w\w)/;
				my $key1 = $1;
				$ncodons{$key1}++;
				$sums{$key1} += $codfreqs{$key};
			}

			foreach my $key (sort keys %codfreqs){
				$key =~ /.{3}(\w\w\w)/;
				my $key1 = $1;
				if ($sums{$key1}==0){
					push @seqfreq, 0;
					push @seqrscu, 0;
				}				
				else {
					push @seqfreq, $codfreqs{$key}/$sums{$key1};
					push @seqrscu, $codfreqs{$key}*$ncodons{$key1}/$sums{$key1};
				}
			}
			$counted = join(' ',@seqfreq);
			$res_freqs{$i} = $id." ".$counted;

			$counted = join(' ',@seqrscu);
			$res_rscus{$i} = $id." ".$counted;
		}

	}

#Print results
    foreach my $key (sort { $a <=> $b } keys %res_counts){
        print COUNT $res_counts{$key}."\n";
    }
	foreach my $key (sort { $a <=> $b } keys %res_freqs){
        print FREQCU $res_freqs{$key}."\n";
    }
	foreach my $key (sort { $a <=> $b } keys %res_rscus){
        print RSCU $res_rscus{$key}."\n";
    }

	return 0;

    close COUNT;
	close FREQCU;
	close RSCU;
	close DATOS;
}

