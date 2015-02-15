#!/usr/bin/perl
use strict;
use warnings;

my $str = " a b c \n ";
chmop($str);
print $str;

sub InesKissMe {
    my ($times,$during) = @_;
    print "Yeah!!! Ines kiss me $times times during $during muites!\n";
}

sub InesKissMeAgain{
    # my $n = 10;
    # @F = [(1) x ($n+1)];
    # print length(@F);
}

$a = 10;
# $d = 2;
# &InesKissMe($a,$d);
# &InesKissMeAgain("Ines Love me!!!!", 1,2,3,4,('go',"fo"),1,222222);

my $seq1 = "abcdef";
print substr($seq1,2,1);

foreach my $z (0..2) {
    print $z,"\n";
}


