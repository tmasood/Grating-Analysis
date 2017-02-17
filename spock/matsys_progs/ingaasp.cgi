#!/usr/bin/perl
print "Content-type:text/html\n\n";

read(STDIN, $buffer, $ENV{'CONTENT_LENGTH'});
@pairs = split(/&/, $buffer);
foreach $pair (@pairs)
{
    ($name, $value) = split(/=/, $pair);
    $value =~ tr/+/ /;
    $value =~ s/%([a-fA-F0-9][a-fA-F0-9])/pack("C", hex($1))/eg;
    $FORM($name) = $value;
}

$lp = $FORM('lp');
    
$c = 2.9979e8;
$h = 6.6261e-34;
$ec = 1.6022e-19;
$ep = ($h*$c)/($lp*1e-6*$ec);

$lambda = $lp;

print "<html><head><title>Form InGaAsP/InP Output</title></head><body>";
print "<h2>Results from InGaAsP/InP Material System </h2> \n"
while ($lambda < 1.8)
{
    # All energies are in Electron Volts
    $a1 = 13.3510 - (5.4554 * $ep) + (1.2332 * ($ep**2));
    $a2 = 0.7140 - (0.3606 * $ep);
	  
    $e1 = 2.5048;
    $e2 = 0.1638;

    # Speed of light 'C', Plank's constant 'H'
    # and electron charge 'EC'

    $e = ($h*$c)/($lambda*1e-6*$ec);
    
    $term1 = $a1/(1 - ((($e/($ep + $e1))**2)));
    $term2 = $a2/(1 - ((($e/($ep + $e2))**2)));
    $epsilon = 1 + $term1 + $term2;
    $effindx = sqrt($epsilon);
    print "$xlambda, $effindx \n";
    $lambda = $lambda + 0.001;
}

print "</body></html>";


