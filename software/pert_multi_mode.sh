#!/bin/sh
#FILE:    pert_multi_mode.sh
#AUTHOR:  Karsten Suhre
#DATE:    Fri Jul 9 14:54:08 CEST 2004
#PURPOSE: perturbate a PDB file following several normal modes
#         compute modes on the fly if not already done
#BUGS:    
#MODIF:   Fri Sep 30 13:31:28 CEST 2005
#         added proj_mode binary file R/W
#         added CUTOFF and NRBL option
#         added diagstd option
#         delete matrix

# elNemo parameters
CUTOFF=8.0
NRBL=0       # set to zero for automatic detection
DIAGSTD="off"

test=${3?"usage: $0 file.pdb [-CUTOFF cutoff] [-NRBL nrbl] [-DIAGSTD]  MODE1 DQ1 MODE2 DQ2 ... "}

# read options (not very elegant ;-)
for i in 1 2 3 ; do
if [ ".$1" = ".-CUTOFF" ] ; then
  CUTOFF=$2
  shift; shift
  echo "using uderdefined CUTOFF=$CUTOFF"
fi
if [ ".$1" = ".-NRBL" ] ; then
  NRBL=$2
  shift; shift
  echo "using uderdefined NRBL=$NRBL"
fi
if [ ".$1" = ".-DIAGSTD" ] ; then
  DIAGSTD="on"
  shift;
  echo "using DIAGSTD"
fi
done


# test again if we still have enough parameters
test=${3?"usage: $0 file.pdb [-CUTOFF cutoff] [-NRBL nrbl]  MODE1 DQ1 MODE2 DQ2 ... "}
fic=$1
command="$*"
shift


if [ ! -f $fic ] ; then
  echo "ERROR: file $fic not found"
  exit 1
fi

base=`basename $fic .pdb`
base=`basename $base .atom`

if [ ! -f $base.eigenfacs ] ; then


  # cleanup first
  rm $base.proj_modes_bin.IO proj_modes_bin.IO pdbmat.sdij? diagrtb.eigenfacs 2> /dev/null

  if [ "$DIAGSTD" = "on" ] ; then # DIAGSTD

    echo "computing normal modes for $base using DIAGSTD"

# run PDBMAT to generate an ASCII file for diagstd
cat >pdbmat.dat <<EOF3
Coordinate FILENAME        = $fic
INTERACtion DISTance CUTOF =      $CUTOFF
INTERACtion FORCE CONStant =      10.000
Origin of MASS values      =       CONS ! CONstant, or from COOrdinate file.
Output PRINTing level      =          2 ! =1: more detailled. =2: debug level.
Bond DEFINITION            =       NONE ! NONe, ALL, or between CONsecutive atoms.
Maximum bond LENGTH        =      0.000
LevelSHIFT                 =    1.0E-09 ! Non-zero value often required (numerical reasons).
Matrix FORMAT              =     FREE   ! Free, or Binary, matrix saved.
EOF3

# generate file matrice.sdijb and matrice.eigenfacs
pdbmat > pdbmat.log || exit 1

# rename matrix to matrice.sdijf for diagstd
mv pdbmat.sdijf matrice.sdijf 

# run diagstd
diagstd > diagstd.log

mv matrice.eigenfacs ${base}.eigenfacs || exit 1
rm matrice.sdijf || exit 1


  else # DIAGRTB

    echo "computing normal modes for $base using DIAGRTB"

  # generate parameter files for pdbmat
  cat > pdbmat.dat <<EOF
 Coordinate FILENAME        = $fic
 INTERACtion DISTance CUTOF =      $CUTOFF
 INTERACtion FORCE CONStant =     10.000
 Origin of MASS values      =       CONS ! CONstant, or from COOrdinate file.
 Output PRINTing level      =          1 ! =1: more detailled. =2: debug level.
 Bond DEFINITION            =       NONE ! NONe, ALL, or between CONsecutive atoms.
 Maximum bond LENGTH        =      0.000
 BOND FORCE CONStant        =      0.000
 ANGLE FORCE CONStant       =      0.000
 LevelSHIFT                 =    1.0E-09 ! Non-zero value often required (numerical reasons).
 Matrix FORMAT              =     BINARY ! FREE ! Free (pdbmat.sdijf), or Binary (pdbmat.sdijb), matrix saved. 
EOF

  # generate matrix
  echo 'running pdbmat'
  pdbmat > pdbmat.log


  # determine NRBL if auto (NRBL=0)
  if [ $NRBL -eq 0 ] ; then
    echo "automatic determination of NRBL (NRBL = nresidues/200 + 1)"
    NONZERO=`grep 'Rdatompdb.*Number of residues found =' pdbmat.log | sed 's/.*found = *//'`
    NRBL=`perl -e "printf \"%i\n\", $NONZERO/200+1;"`
    echo "$NONZERO non-zero elements, NRBL set to $NRBL"
  fi

  # generate parameter files for diagrtb
  cat > diagrtb.dat <<EOF2
 MATRix filename            = pdbmat.sdijb ! 
 COORdinates filename       = $fic
 Eigenvector OUTPut filename= diagrtb.eigenfacs
 Nb of VECTors required     =        106
 EigeNVALues chosen         =       LOWE   ! LOWEst, HIGHest.
 Type of SUBStructuring     =       NONE   ! RESIdues, SECOndary, SUBUnits, DOMAins, NONE.
 Nb of residues per BLOck   =          $NRBL
 Origin of MASS values      =       CONS   ! CONStant, COORdinate, PDB.
 Temporary files cleaning   =       ALL    ! ALL, NOne.
 Matrix FORMAT              =     BINARY ! FREE ! Free (pdbmat.sdijf), or Binary (pdbmat.sdijb), matrix saved. 
 Output PRINting level      =          1   ! =1: More detailed; =2: Debug level.
EOF2

  # get eigenvectors
  echo 'running diagrtb'
  diagrtb > diagrtb.log

  # rename files
  rm pdbmat.sdijb  || exit
  mv diagrtb.eigenfacs $base.eigenfacs

  fi # DIAGRTB

fi

# check if all went right
if [ ! -f $base.eigenfacs ] ; then
  echo "ERROR: normal mode computation failed!"
  exit 1
fi


# perturbate with normal modes
tmp=eraseme.$$.pdb
cp $fic $tmp
while [ $# -ge 2 ] ; do
  MODE=$1
  DQ=$2
  echo "applying MODE=$MODE with DQ=$DQ to $base"
  shift ; shift

  # run proj_modes
  rm projmod_dq.pdb 2> /dev/null
  # use binary mode
  # proj_modes > proj_modes.log <<EOF
  if [ -e $base.proj_modes_bin.IO ] ; then
    echo "using binary file $base.proj_modes_bin.IO"
    mv $base.proj_modes_bin.IO proj_modes_bin.IO
  fi
  proj_modes_bin > proj_modes.log  <<EOF
$base.eigenfacs
$tmp
nimportequoi
n
y
$DQ
$MODE
EOF
  mv proj_modes_bin.IO $base.proj_modes_bin.IO
  mv projmod_dq.pdb $tmp

done

# generate output file
echo "REMARK NMA PERTURBED MODEL" >$base.NMApert.pdb
echo "REMARK COMMAND LINE: \"$0 $command\"" >>$base.NMApert.pdb
cat $tmp | perl -e '
#FILE:    pdb_fix.pl
#AUTHOR:  Karsten Suhre
#DATE:    Thu Jun 10 10:16:33 CEST 2004
#PURPOSE: set occupancy to one
#         set B-factors to 20.00
#         renumber ATOMs
#         add chain ids (every time when resnum decreases
#         
#BUGS:    
#MODIF:   


# "Beware of bugs in the above code; I have only proved it correct, not
# tried it."
#               -- Donald Knuth

$resnumold = -9999;
@CHAIN = ("A".."Z", "0", "1".."9");
$nchain = 0;
$nmulti = 0;
$natom = 1;
while (<>) {
  if ( ( /^HETATM/ ) or ( /^ATOM/ ) ) {

    # treat chain ids
    $resnum = substr($_,22,4);
    if ($resnum < $resnumold) { # start a new chain
      $nchain++;
      if ($nchain>$#CHAIN) {$nchain = 0; $nmulti++};
    }
    $resnumold = $resnum;

    printf substr($_,0,6);
    printf "%5i", $natom;
    print substr($_,11,10);
    print $CHAIN[$nchain];
    print substr($_,22,31);
    print "  1.00 20.00\n";
    $natom++;
  } else {
    print;
  }
}
' >>$base.NMApert.pdb
rm $tmp
# cleanup
rm pdbmat.xyzm pdbmat.dat_run pdbmat.dat diagrtb.dat_run diagrtb.dat diagrtb.eigenfacs 2> /dev/null
exit 0
