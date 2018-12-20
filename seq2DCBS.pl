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
my $file1;
my %res_DCBS;



if (@ARGV){
	$file1 = $ARGV[0];
	} 
else{
    print "Please type in the file name: ";
    $file1 = <STDIN>;
}


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
	$head = $head." ".$_
}


$file1 =~ /(^.+)\./;
my $name = $1;
print $name;
my $file2 = $name.".DCBS";

#input.dat and temp.file
$tempfile = readfile ($file1,"seq.temp");

#calculate DCBS
DCBS ($tempfile,$file2);

unlink "seq.temp";
unlink "freq.tmp";


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
sub DCBS {
	#$tempfile,
	my @files=@_;
	my $file=$files[0];
	my $dcbs=$files[1];
	open (TEMP,"<$file");
	open (RESULTADOS,">$dcbs");

    $res_DCBS{0} = $head." DCBS";

	chomp(my @lines = <TEMP>);

	my $id;
	my $i=0;

	foreach (@lines){
		my $line = $_;
		my %codoncount;
		my %codoncount2;
		my @seqDCBS;
		my %basesfreq;
		my %basesposcount;
		my %basesposfreq;
		my %codfreqs;
		my %sums;
		my @codonorder;
		my %fxyz;
		my %dxyz;

		if ($line =~ /^>/){
			$line =~ s/>//;
			$id = $line;
	        $i++;
		}
		else {
			my $baselen = length $line;
            my $codonlen = ($baselen-6)/3;
			@seqcodons = $line =~ /(...)/g;
            shift @seqcodons;
            pop @seqcodons;
			my $bpos;
			# calculate f(x) f(y) f(z) gloabl and possition speciffic
			foreach my $codon (@seqcodons){
			  	my @xyz = $codon =~ /(.)/g;
				for (my $i=0; $i < 3; $i++) {
					$bpos = $xyz[$i].$i;
					$basesposcount{$bpos}++;
			 	}
			}
                       
			foreach my $key (sort keys %basesposcount){
                $basesposfreq{$key}=0;
				$basesposfreq{$key}=$basesposcount{$key}/$codonlen; 
    		}

			# calcualte codon counts
			my %codoncount;
			my %codoncount2;
			foreach my $acod (sort keys %cod){
				$codoncount{$acod} = 0;
				$codoncount2{$cod{$acod}}=0;
			}

			foreach my $codon (@seqcodons){
				$codoncount2{$codon}++;
			}
			foreach (sort keys %cod){
				my $acod=$_;
				$codoncount{$acod} = $codoncount2{$cod{$acod}};
			}

			#calculate freqs
			foreach my $key (sort keys %cod){
				push @codonorder,$key;
				}

			for (my $i=0; $i <= 63; $i++){
				$codfreqs{$codonorder[$i]}=$codoncount{$codonorder[$i]};

			}

			foreach my $key (sort keys %codfreqs){
				$fxyz{$key} = $codfreqs{$key}/$codonlen;
				}
			

			#Calculate dxyz, for fxyz = 0, different approaches were made. Firts, the dxyz value
            #was set to 0
			foreach my $key (sort keys %fxyz){
				$key =~ /.{7}(\w\w\w)/;
				my $key1 = $1;
				my @pos = $key1 =~ /(.)(.)(.)/g;
                my $bigger;
                my $smaller;
                if ($fxyz{$key} != 0){
				    $bigger = $fxyz{$key}/($basesposfreq{$pos[0]."0"}*$basesposfreq{$pos[1]."1"}*$basesposfreq{$pos[2]."2"});
				    $smaller = ($basesposfreq{$pos[0]."0"}*$basesposfreq{$pos[1]."1"}*$basesposfreq{$pos[2]."2"})/$fxyz{$key};
                    
				    if ($bigger > $smaller){
					    $dxyz{$key} = $bigger;
				    }
				    else {
					    $dxyz{$key} = $smaller;
				    }
                }
                else{
                    $dxyz{$key} = 0;
                }
                
			}
            #my @higher = sort { $dxyz{$a} <=> $dxyz{$b} } keys %dxyz;
            # for codons with fxyz = 0 != Met, trp, TER, assigns to max(dxyz)+10 -- Anyway it
            # doesn't affect DCBS 
            #foreach my $key (sort keys %dxyz){
            #    if ($dxyz{$key} == 0){
            #       $dxyz{$key} = $dxyz{$higher[-1]}+10;
            #    }
            #}

            foreach my $key (sort keys %dxyz){
                if ($key =~ /TER/){
                   $dxyz{$key} = 0;
                }
                elsif ($key =~ /Trp/){
                   $dxyz{$key} = 0;
                }
                elsif ($key =~ /Met/){
                   $dxyz{$key} = 0;
                }

 		        push @seqDCBS, $dxyz{$key};
            }
            #calculates DCBS
            my $DCBS;
            foreach my $codon (@seqcodons){
                foreach my $key (sort keys %dxyz){
                    if ($key =~ $codon){
                    $DCBS += $dxyz{$key};
                    }
                }
            }
            $DCBS = $DCBS/$codonlen;
			push @seqDCBS, $DCBS;
			
			my $counted = join(' ',@seqDCBS);
			$res_DCBS{$i} = $id." ".$counted;
		}
	}

	foreach my $key (sort { $a <=> $b } keys %res_DCBS){
        print RESULTADOS $res_DCBS{$key}."\n";
    }

	return;
	close TEMP;
	close RESULTADOS;
}

