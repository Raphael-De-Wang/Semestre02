#!/usr/bin/perl

# Needleman-Wunsch algorithm
# ----------------
# Inspired from:
# Peter Sestoft, Royal Veterinary and Agricultural University, Denmark
# Reference: http://www.dina.kvl.dk/~sestoft/bsa.html
# sestoft@dina.kvl.dk * 2003-04-19, 2003-05-04, 2003-08-25, 2003-10-16
# ----------------

use strict;
use warnings;
use myConstants;
use CGI qw/:standard/;  # to use param()

# ----------------------------------------------------------------------
# Set global variables

my $true = 1;
my $false= 0;

# $seq1, $seq2: The sequences to align (default sequences given)
my ($seq1, $seq2) = ("TACGT","CTAGGAT");

# Set nucleotide substitution matrix
my $matrix = \@Constants::matchMismatch;

# Set the gap extension value
# my $e = 2;
my $wg= 5;
my $ws= 1;

print ("\n# My sequences:\n# SEQ1 --> ".$seq1."\n# SEQ2 --> ".$seq2."\n"); 

# ----------------------------------------------------------------------
# Global constants
my $sepString = "# --------------------------------\n";
print($sepString);

# The traceback matrices are indexed by (direction, row, column).
my @DIR = (1, 2, 3);

# Directions in the linear (2D) traceback: 
# 0=stop; 1=from North (above); 2=from Northwest; 3=from West (left)
my ($FROMN, $FROMNW, $FROMW) = @DIR;

# Color codes for the traceback
my ($RED, $BLUE, $GREEN) = (1, 2, 3);


# ----------------------------------------------------------------------
# The Needleman-Wunsch global alignment algorithm
# Input: The sequences $x and $y to align, the amino acids sub matrix
# Output: Output: references to the F and B matrices, and the aligned sequences
sub globalAlignLinear {
    # Get the parameters (By ref, val, val)
    my ($matrix, $x, $y) = @_;
    
    # Get the length of the sequences
    my ($n, $m) = (length($x), length($y));
    
    # The dynamic programming matrix
    my @V;
    for (my $j=0; $j<=$m; $j++) {
        $V[$j] = [(0) x ($n+1)];
    }
    
    # The traceback matrix
    my @B;
    foreach my $dir (@DIR) {
        for (my $j=0; $j<=$m; $j++) {
            $B[$dir][$j] = [(0) x ($n+1)];
        }
    }

    # S ------i-------
    # T -------------j 
    my @E;
    for (my $j=0; $j<=$m; $j++) {
        $E[$j] = [(0) x ($n+1)];
    }

    # S --------------i 
    # T ------j--------
    my @F;
    for (my $j=0; $j<=$m; $j++) {
        $F[$j] = [(0) x ($n+1)];
    }
    
    # S ------i 
    # T ------j 
    my @G;
    for (my $j=0; $j<=$m; $j++) {
        $G[$j] = [(0) x ($n+1)];
    }

    # Initialize upper and left-hand borders of V and B matrices
    for (my $i=1; $i<=$n; $i++) {
	$E[0][$i] = -$wg - $ws * $i;
        $V[0][$i] = -$wg - $ws * $i;
        $B[$FROMW][0][$i] = $RED;
    }
    
    for (my $j=1; $j<=$m; $j++) {
	$F[$j][0] = -$wg - $ws * $j;
        $V[$j][0] = -$wg - $ws * $j;
        $B[$FROMN][$j][0] = $RED;
    }

    # Do the iteration
#######################################
    my ($maxScore,$tN,$tNW,$tW) = (0,0,0,0);
    for (my $i=1; $i<=$n; $i++) {
	for (my $j=1;$j<=$m; $j++) {
	    ($maxScore,$tN,$tNW,$tW) = &bestScore($matrix,$wg,$ws,\@V,\@E,\@F,\@G,$x,$y,$i,$j);
	    $V[$j][$i] = $maxScore;
	    
	    $B[$FROMN][$j][$i]  = ( $tN  == $false ? 0 : 1 ) ;
	    $B[$FROMNW][$j][$i] = ( $tNW == $false ? 0 : 1 ) ;
	    $B[$FROMW][$j][$i]  = ( $tW  == $false ? 0 : 1 ) ;
	}
    }
#######################################
    &markReachable2(\@B, $n, $m);
    return (\@V, \@B, &traceback2($x, $y, \@B, $n, $m));  
}

# ----------------------------------------------------------------------
# Common subroutines for linear gap cost routines

# Reconstruct the alignment from the traceback, backwards, from ($i, $j)
sub traceback2 {
  my ($x, $y, $B, $i, $j) = @_;         # B by reference
  my ($xAlign, $yAlign) = ("", "");
  while ($$B[$FROMN][$j][$i] || $$B[$FROMW][$j][$i] || $$B[$FROMNW][$j][$i]) {
    if ($$B[$FROMN][$j][$i]) {
      $$B[$FROMN][$j][$i] = $GREEN;
      $xAlign .= "-"; 
      $yAlign .= substr($y, $j-1, 1);
      $j--;
    } elsif ($$B[$FROMW][$j][$i]) {
      $$B[$FROMW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= "-"; 
      $i--;
    } elsif ($$B[$FROMNW][$j][$i]) {
      $$B[$FROMNW][$j][$i] = $GREEN;
      $xAlign .= substr($x, $i-1, 1);
      $yAlign .= substr($y, $j-1, 1);
      $i--; $j--;
    }
  }
  # Warning: these expressions cannot be inlined in the list
  $xAlign = reverse $xAlign;
  $yAlign = reverse $yAlign;
  return ($xAlign, $yAlign);
}

# Mark all traceback arrows reachable from a ($i, $j)
sub markReachable2 {
  my ($B, $i, $j) = @_;         # B by reference
  if ($$B[$FROMN][$j][$i] == $RED) {
    $$B[$FROMN][$j][$i] = $BLUE;
    &markReachable2($B, $i, $j-1);
  } 
  if ($$B[$FROMW][$j][$i] == $RED) {
    $$B[$FROMW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j);
  } 
  if ($$B[$FROMNW][$j][$i] == $RED) {
    $$B[$FROMNW][$j][$i] = $BLUE;
    &markReachable2($B, $i-1, $j-1);
  }
}

# ----------------------------------------------------------------------
# Main program: 
# Check arguments and call the alignment and graph drawing functions

my $num_args = $#ARGV + 1;
if ($num_args >= 2) {
    $seq1 = $ARGV[0];
    $seq1 = $ARGV[1];
}

print( "## -- Global Alignement -- ##\n" );
my ($F, $B, $xa, $ya) = &globalAlignLinear($matrix, $seq1, $seq2);

print("--> $xa\n--> $ya\n\n# Dynamic Programming Matrix\n");

my $cellSize = 5; 
for (my $i=0; $i<=length($seq2); $i++) {
    for (my $j=0; $j<=length($seq1); $j++) {
        for( my $k = length($$F[$i][$j]); $k <= $cellSize; $k++) {
            print(" ");
        }
        print( $$F[$i][$j] ); 
    }
    print("\n");
}

exit();

print("\n\n Arrows color code:
            --> RED(1): Possible
            --> BLUE(2): Reachable in traceback
            --> GREEN(3): Choosen in traceback");
print("\n\n# Arrows North\n");
$cellSize = 5; 
for (my $i=0; $i<=length($seq2); $i++) {
    for (my $j=0; $j<=length($seq1); $j++) {
        for( my $k = length($$B[$FROMN][$i][$j]); $k <= $cellSize; $k++) {
            print(" ");
        }
        print( $$B[$FROMN][$i][$j] ); 
    }
    print("\n");
}

print("\n\n# Arrows West\n");
$cellSize = 5; 
for (my $i=0; $i<=length($seq2); $i++) {
    for (my $j=0; $j<=length($seq1); $j++) {
        for( my $k = length($$B[$FROMW][$i][$j]); $k <= $cellSize; $k++) {
            print(" ");
        }
        print( $$B[$FROMW][$i][$j] ); 
    }
    print("\n");
}

print("\n\n# Arrows North West\n");
$cellSize = 5; 
for (my $i=0; $i<=length($seq2); $i++) {
    for (my $j=0; $j<=length($seq1); $j++) {
        for( my $k = length($$B[$FROMNW][$i][$j]); $k <= $cellSize; $k++) {
            print(" ");
        }
        print( $$B[$FROMNW][$i][$j] ); 
    }
    print("\n");
}
