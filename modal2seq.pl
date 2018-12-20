#! /usr/bin/env perl

use strict;
use warnings;

#declaración de variables
my $freq;
my $name1;
my $codons;
my @codon_a;
my %cod;

# abre archivos de codons, <FH> en contexto escalar lee una linea. si se pone otra vez -> la siguiente
open (CODONS, "<codons.txt");
$codons = <CODONS>;
close (CODONS);
@codon_a = split(",", $codons);

if (@ARGV){
	$freq = $ARGV[0];
	open(FREQ, "<$freq");
	$freq =~ /(^.+)\./;
	$freq=$1;
	$name1 = $ARGV[1];
	open (RESULTADOS,">$name1");
} 
else{
	print "Introduzca el nombre del archivo de frecuencias modales\n";
	$freq = <STDIN>;
	open(FREQ, "<$freq");
	$freq =~ /(^.+)\./;
	$freq=$1;
	$name1 = $freq.".seq";	
	open (RESULTADOS,">$name1");
}

#declarar variables de todo el namespace-package


while (<FREQ>) {
	my $line = $_;
	chomp $line;
	$line =~ s/^>//;
	$line =~ s/\|/,/g;
	(my $lineid, my $modal) = split("\t",$line);
	my @modal_a = split(",",$modal);
	my $naa = 100;
	foreach my $x (@modal_a) { 
		$x = $x * $naa; 
	}	
	for (my $i=0; $i <= 58; $i++){
		$cod{$codon_a[$i]} = $modal_a[$i];
	}
	#algoritmo
	print RESULTADOS ">".$freq."\n"."ATG";
	foreach my $key (keys %cod){
	print RESULTADOS "$key"x$cod{$key};
	}
	print RESULTADOS "TGA\n";
}
close FREQ;
close RESULTADOS;

