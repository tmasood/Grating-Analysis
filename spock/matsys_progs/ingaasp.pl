#!/usr/bin/perl

$xlambda = 0.93;
while ($xlambda < 1.8)
{
    # All energies are in Electron Volts
    $ep = 1.32;
    $a1 = 13.3510 - (5.4554 * $ep) + (1.2332 * ($ep**2));
    $a2 = 0.7140 - (0.3606 * $ep);
	  
    $e1 = 2.5048;
    $e2 = 0.1638;

    # Speed of light 'C', Plank's constant 'H'
    # and electron charge 'EC'
    $c = 2.9979e8;
    $h = 6.6261e-34;
    $ec = 1.6022e-19;
    $e = ($h*$c)/($xlambda*1e-6*$ec);
    
    $term1 = $a1/(1 - ((($e/($ep + $e1))**2)));
    $term2 = $a2/(1 - ((($e/($ep + $e2))**2)));
    $epsilon = 1 + $term1 + $term2;
    $effindx = sqrt($epsilon);
    print STDOUT "$xlambda, $effindx \n";
    $xlambda = $xlambda + 0.001;
}


