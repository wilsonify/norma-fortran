#!/bin/sh
#FILE:    gen_pert.sh
#AUTHOR:  Karsten Suhre
#DATE:    Thu Jun 23 14:07:07 CEST 2005
#PURPOSE: generate a series of NMA pertured models 
#         the resulting anim.pdb may be best viewed using vmd
#BUGS:    
#MODIF:   

line="${1?'usage: gen_pert.sh mode1 dq1 mode2 dq2 ... modeN dqN'}"

> anim.pdb
for i in 0 1 2 3 4 5 6 7 8 9 10 ; do

  # generate perturbed model
  echo $i " $*" | perl -e '
$_ = <>;
@d = split /  */;
for ($i=1; $i<$#d; $i+=2) {
  print $d[$i], " ", $d[$i+1]*0.1*$d[0], " ";
}
print "\n";
' > func.in

  cat func.in
  sh func.sh
  grep '#' `ls -tr RUN.minimize/o/f*.s | tail -1`

  echo "MODEL $i" >> anim.pdb
  grep ATOM RUN.minimize/Aa.pdb >> anim.pdb
  echo "ENDMDL $i">> anim.pdb

done

echo "Your animated PDB file is anim.pdb"
