# ----------------
# Inspired from:
# Peter Sestoft, Royal Veterinary and Agricultural University, Denmark
# Reference: http://www.dina.kvl.dk/~sestoft/bsa.html
# sestoft@dina.kvl.dk * 2003-04-19, 2003-05-04, 2003-08-25, 2003-10-16
# ----------------

use strict;
use warnings;
use List::Util qw/sum/;

# The amino acids, and a hashmap from amino acid to its index 0-19
@Constants::aas  = split(//, "ACTG");
my $len = 4;

my %aaIndex;
for (my $i=0; $i<$len; $i++) {
  $aaIndex{$Constants::aas[$i]} = $aaIndex{lc($Constants::aas[$i])} = $i;
}

# The score of a pair of amino acids in a given matrix
sub score {
  my ($matrix, $aa1, $aa2) = @_;	# By ref, val, val
  return $$matrix[$aaIndex{$aa1}][$aaIndex{$aa2}];
}

# The maximum of a list of numbers
sub max {
  my $res = -1E300;
  foreach (@_) {
    if ($_ > $res) {
      $res = $_;
    }
  }
  return $res;
}

# ----------------------------------------------------------------------
# Calcule the best score for F[i,j]
# Input: ScoreMatrix, gap extension value, The dynamic programming matrix, 2 sequences, position i,j
# Output: best score for F[i,j]
my $true = 1;
my $false= 0;
sub bestScore {
    my ($matrix,$wg,$ws,$V,$E,$F,$G,$x,$y,$i,$j) = @_;
    my $scoreNW= $$V[$j-1][$i-1] + score($matrix,substr($x,$i-1,1),substr($y,$j-1,1)); # match/mismatch
    my $scoreN = $$V[$j-1][$i] - $wg - $ws; # add a gap
    my $scoreW = $$V[$j][$i-1] - $wg - $ws; # add a gap

    $$G[$j][$i] = $scoreNW;
    $$E[$j][$i] = max(($$E[$j-1][$i]-$ws,$scoreN));
    $$F[$j][$i] = max(($$F[$j][$i-1]-$ws,$scoreW));
    
    # my $maxScore = max(($scoreN,$scoreNW,$scoreW));
    my $maxScore = max(($$E[$j][$i],$$F[$j][$i],$$G[$j][$i]));
    
    my $tNW= $maxScore==$$G[$j][$i]?$true:$false;
    my $tN = $maxScore==$$E[$j][$i]?$true:$false;
    my $tW = $maxScore==$$F[$j][$i]?$true:$false;
    
    return ($maxScore,$tN,$tNW,$tW);
}

# A simple match/mismatch nucleotide substitution matrix
@Constants::matchMismatch = 
    #  A  T  C  G  
  ( [  2,-1,-1,-1],  # A
    [ -1, 2,-1,-1],  # T
    [ -1,-1, 2,-1],  # C
    [ -1,-1,-1, 2]   # G
    #  A  T  C  G  
    );

