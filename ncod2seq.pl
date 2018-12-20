#! /usr/bin/env perl

use strict;
use warnings;

#declaración de variables
my @freq_a;
my $freq;
my $name1;

# abre archivos
if (@ARGV){
	$freq = $ARGV[0];
	open(FREQ, "<$freq");
	@freq_a = <FREQ>;
	close FREQ;
	$name1 = $ARGV[1];
	open (RESULTADOS,">$name1");
} 
else{
	print "Introduzca el nombre del archivo de frecuencias\n";
	$freq = <STDIN>;
	open(FREQ, "<$freq");
	@freq_a = <FREQ>;
	close FREQ;
	print "Introduzca el nombre del archivo de resultados\n";
	$name1 = <STDIN>;
	open (RESULTADOS,">$name1");
}

#declarar variables de todo el namespace-package
my %freq_h;

#map to hash split(",", $_);
%freq_h = map { split( "\t", $_) } @freq_a;
my @seqname = split /\./, $name1;

#algoritmo
print RESULTADOS ">".$seqname[0]."\n"."ATG";
foreach my $key (keys %freq_h){
print RESULTADOS "$key"x$freq_h{$key};
print RESULTADOS "TGA";
}
close RESULTADOS;

