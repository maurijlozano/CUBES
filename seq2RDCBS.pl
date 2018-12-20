#! /usr/bin/env perl
#Programed by Mauricio J Lozano
#La idea es calcular DCBS relativo al DCBS de la modal de las PHE
#en memoria en vez de leyendo y excribiendo, pero es mas lento...
use strict;
use warnings;
$| = 1;

#variables
my %cod;
my %aa;
my @codons_a;
my @seqcodons;
my @tempfile;
my @mtempfile;
my $counts;
my $head;
my $file1;
my $modal;
my @resultados;

if (@ARGV){
	$file1 = $ARGV[0];
    $modal = $ARGV[1];
	} 
else{
    print "Please type in the file name: ";
    $file1 = <STDIN>;
    print "Please type in the modal file name: ";
    $modal = <STDIN>;
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
	$head = $head." ".$_;
}


#input.dat and temp.file
	
@tempfile = readfile ($file1);
@mtempfile = readfile ($modal);

my $temp ="seq.temp";
open (TEMP,">$temp");
print TEMP @tempfile;
close TEMP;

my $mtemp ="modal.temp";
open (TEMP,">$mtemp");
print TEMP @mtempfile;
close TEMP;

$file1 =~ /(^.+)\./;
my $name = $1;
print $name;
my $file2 = $name.".RDCBS";

#calculate DCBS
@resultados = RDCBS ($temp,$mtemp);



open (RESULTADOS,">$file2");
print RESULTADOS @resultados;
close RESULTADOS;


unlink "seq.temp";
unlink "modal.temp";
#**********************************************
sub readfile {
	my $input=$_[0];
    my @temp;
	open (DATOS, "<$input");
    chomp(my @lines = <DATOS>);
    foreach my $line (@lines) {
		if ($line =~ /^>/){
			$line =~ s/\s/_/g;
			$line =~ s/\'/_/g;
			push @temp, "\n".$line."\n";
		}
		elsif($line =~ /^$/){
		}
		else {
			$line =~ s/\r//g;
            $line =~ s/-*$//;
			$line =~ tr/agtc/AGTC/;
			push @temp, $line;
			}
	}
	close DATOS;
	return @temp;
}

sub RDCBS {
   	my $file=$_[0];
    my $mfile=$_[1];
	open (TEMP,"<$file");
	open (MTEMP,"<$mfile");
    chomp(my @lines = <TEMP>);
    chomp(my @mlines = <MTEMP>);
    shift @lines;
    shift @mlines;
    my @results;
    push @results, $head." RDCBS\n";
	
    my %mdxyz;
    foreach my $line (@mlines) {
		my %codoncount;
		my %basesfreq;
		my %basesposcount;
		my %basesposfreq;
		my %codfreqs;
		my %sums;
		my @codonorder;
		my %fxyz;
		my %dxyz;
		
        chomp $line;
		if ($line =~ /^>/){
			}
        elsif ($line =~ /^\n/){
		    }
		else {
            chomp $line;
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
				$basesposfreq{$key}=$basesposcount{$key}/$codonlen;
			}

			# calcualte codon counts
			foreach my $codon (@seqcodons){
				foreach (sort keys %cod){
					my $acod=$_;
					if ($codon =~ /^$cod{$acod}$/i){
						$codoncount{$acod}++;
					}
					if (not exists $codoncount{$acod}) {
						$codoncount{$acod}=0;
					}
				}
			}


			foreach my $key (sort keys %cod){
				push @codonorder,$key;
				}
			#calculate freqs
			for (my $i=0; $i <= 63; $i++){
				$codfreqs{$codonorder[$i]}=$codoncount{$codonorder[$i]};
			}
			foreach my $key (sort keys %codfreqs){
			
					$fxyz{$key} = $codfreqs{$key}/$codonlen;
				}
			
			#Calculate dxyz, for fxyz = 0, different approaches were made. Firts, the dxyz value was set to 0
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
           #calcula la media aritmética de los valores != 0
           my $dxyzam;
           my $count;
           foreach my $key (sort keys %dxyz){
                    if ($dxyz{$key} != 0) {
                       $dxyzam += $dxyz{$key};
                       $count++;
                    }
                }
           $dxyzam=$dxyzam/$count;

           foreach my $key (sort keys %dxyz){
 		       if ($dxyz{$key} == 0 ){
                    $mdxyz{$key} = $dxyzam;
                }
                else {
                     $mdxyz{$key} = $dxyz{$key};
                }

            }
        }
    }



    foreach my $line (@lines) {
		my %codoncount;
		my %basesfreq;
		my %basesposcount;
		my %basesposfreq;
		my %codfreqs;
		my %sums;
		my @codonorder;
		my %fxyz;
		my %dxyz;
	
        chomp $line;
		if ($line =~ /^>/){
			$line =~ s/>//;
			push @results, $line." ";
		}
        elsif ($line =~ /^\n/){
		    }
		else {
            chomp $line;
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
			foreach my $codon (@seqcodons){
				foreach (sort keys %cod){
					my $acod=$_;
					if ($codon =~ /^$cod{$acod}$/i){
						$codoncount{$acod}++;
					}
					if (not exists $codoncount{$acod}) {
						$codoncount{$acod}=0;
					}
				}
			}


			foreach my $key (sort keys %cod){
				push @codonorder,$key;
				}
			#calculate freqs
			for (my $i=0; $i <= 63; $i++){
				$codfreqs{$codonorder[$i]}=$codoncount{$codonorder[$i]};
			}
			foreach my $key (sort keys %codfreqs){
				$fxyz{$key} = $codfreqs{$key}/$codonlen;
				}

			#Calculate dxyz, for fxyz = 0, different approaches were made. Firts, the dxyz value was set to 0
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

 		        push @results, $dxyz{$key}/$mdxyz{$key}." ";
                $dxyz{$key} = $dxyz{$key}/$mdxyz{$key};
            }
#calculates DCBS
            my $RDCBS;
            foreach my $codon (@seqcodons){
                foreach my $key (sort keys %dxyz){
                    if ($key =~ $codon){
                    $RDCBS += $dxyz{$key};
                    }
                }
            }
            $RDCBS = $RDCBS/$codonlen;
            push @results, $RDCBS."\n";
		}
	}
	return @results;
}

