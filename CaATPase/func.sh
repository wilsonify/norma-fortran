#!/bin/sh
#FILE:    func.sh
#AUTHOR:  Karsten Suhre
#DATE:    Fri Jul 9 15:54:25 CEST 2004
#PURPOSE: run NMA and URO using input func.in from NORMA.exe and output func.out
#BUGS:    
#MODIF:   


# "Beware of bugs in the above code; I have only proved it correct, not
# tried it."
#               -- Donald Knuth

# redirect all output to log file
exec > func.log 2>&1

# the name of the PDB file
pdb=1eul.moved


base=`basename $pdb .pdb`

# clean-up old files 
rm func.out $base.NMApert.pdb 2> /dev/null

# do some checks
if [ ! -f func.in ] ; then
  echo "ERROR: no input parameters: func.in"
  exit 1
fi

if [ ! -e ./d/$base.pdb ] ; then
  echo "ERROR: no PDB file: ./d/$base.pdb"
  exit 1
fi

echo "func.sh computing NMA and URO for input parameters "`head -1 func.in`

# generate a perturbed model using NMA
echo pert_multi_mode.sh ./d/$base.pdb `head -1 func.in`
pert_multi_mode.sh ./d/$base.pdb `head -1 func.in`

cat > ./d/$base.NMApert.pdb <<EOF
FORMAT (12X,A4,14X,3F8.3,6X,F6.2)                                       
EOF
cat $base.NMApert.pdb | grep ATOM >> ./d/$base.NMApert.pdb
rm $base.NMApert.pdb

# run URO
./go.URO.minimize || exit 1

echo "func.sh terminated for input parameters "`head -1 func.in`
