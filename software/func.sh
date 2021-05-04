#!/bin/sh
#FILE:    func.sh
#AUTHOR:  Karsten Suhre
#DATE:    Fri Sep 23 12:12:04 CEST 2005
#PURPOSE: a dummy func.sh to test NORMA.f
#         see shell code below (after "exit" for a more realistic func.sh)
#BUGS:    
#MODIF:   


exec > func.log 2>&1

rm func.out 2> /dev/null

# a dummy function to minimize (replace this by a call to NMA and URO - for an example see below)
head -1 func.in | perl -e '
$in = <>;
chomp $in; $in =~ s///; $in =~ s/^ *//;;
@d = split /  */, $in;
print "IN: $in\n";
print $#d + 1, " elements \n";
$f = 0.0;
$x = 0.0;
for ($i=0; $i<$#d; $i+=2) {
  print "$d[$i] $d[$i+1]\n";
  $f += $d[$i+1] * $d[$i+1];
  $x += $d[$i+1];
}
$f += 100*(1-cos ($x/1.));
$f *= 10;
print "f=$f\n";
system "echo $f > func.out";
'

cat func.out

