#FILE:    test.sh
#AUTHOR:  Karsten Suhre
#DATE:    Fri Apr 7 12:18:48 CEST 2006
#PURPOSE: test installation of NORMA
#BUGS:    
#MODIF:   

# test option (all = run all tests , this is the default)
# individual test options: NORMA, groel, ibdv, CaATPase
opt=${1-all}

echo; echo; echo
echo " ------ NORMA INSTALLATION TEST STARTED ------ "
echo
echo "This test will perform a number of computations using the different parts of NORMA."
echo "Testing component: $opt"
echo "It will abort if an unexpected error occurs."
echo "Should you be unable to solve the problem, please contact"
echo "Karsten Suhre <karsten.suhre@igs.cnrs-mrs.fr>"
echo

# running setup script
if [ ! -f ./software/setup.sh ] ; then

  echo "ERROR: NORMA is not properly installed (missing file ./software/setup.sh)"
  exit 1

fi
#source ./software/setup.sh
./software/setup.sh

#------------------------------------------------------------------------
# test NORMA
#------------------------------------------------------------------------
if [ "$opt" == "all" -o "$opt" == "NORMA" ] ; then
cd ./software || exit 1
echo "testing NORMA minimization code (amebsa) .."
NORMA.exe  > NORMA.log
echo "testing NORMA minimization code (amebsa) .. done"

echo "testing normal mode code .. "
pert_multi_mode.sh 1hhp.pdb 11 10.0 12 -15.0 > pert_multi_mode.log
echo "testing normal mode code .. done"

echo
echo "The following two lines should correspond :"
echo " Final minimum =   <a small number; about 1E-4>"
tail -1 NORMA.log

echo
echo "The following two lines should correspond (when compiling with g77) :"
echo "ATOM      1  N   PRO A   1      52.688  59.090  -8.10  1.00 20.00"
head -3 1hhp.NMApert.pdb | tail -1

echo "checking URO installation .."
if [ -e $URO/setup ] ; then
  echo "URO is installed in directory $URO"
else
  cat <<'EOF1'
  ERROR: URO not found make sure that URO is installed and that the environment
         variable $URO points to URO installation directory.
EOF1
  exit 1
fi
echo "checking URO installation .. done"
cd ..
fi
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# test groel case
#------------------------------------------------------------------------
if [ "$opt" == "all" -o "$opt" == "groel" ] ; then
echo "testing groel case .."
cd groel || exit 1

# setup URO in groel directory
echo "setting up URO in groel directory"
csh $URO/setup || exit 1
mv sym gs.sym gs.real symlist GroEl.ezd 1aonA.pdb ./d 

# run a simple function evaluation
echo "running one step of NMA and URO on groel"
sh func.sh 
cat func.log

# run a one-mode minimisation
echo "Full fitting with NORMA (based on NORMA.inp)"
time NORMA.exe || exit 1

# generate an animated pdb file for the solution
echo "computing an animation for the final NORMA fitting"
mv func.in final.in
gen_pert.sh `cat final.in`

echo "testing groel case .. done"
cd ..
fi
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# test ibdv case
#------------------------------------------------------------------------
if [ "$opt" == "all" -o "$opt" == "ibdv" ] ; then
echo "testing ibdv case .."
cd ibdv || exit 1

# setup URO in ibdv_vp2 directory
echo "setting up URO in ibdv directory"
csh $URO/setup || exit 1
mv sym gs.sym gs.real symlist map_odd.ezd 1WCD.reo.pdb ./d 

# run a simple function evaluation
echo "running one step of NMA and URO on ibdv"
sh func.sh 
cat func.log

echo "testing ibdv case .. done"
cd ..
fi
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# test CaATPase case
#------------------------------------------------------------------------
if [ "$opt" == "all" -o "$opt" == "CaATPase" ] ; then
echo "testing CaATPase case .."
cd CaATPase || exit 1

# setup URO in CaATPase directory
echo "setting up URO in CaATPase directory"
csh $URO/setup || exit 1
mv sym gs.sym gs.real symlist dtg8.ezd 1eul.moved.pdb ./d 

# run a simple function evaluation
echo "running one step of NMA and URO on CaATPase"
sh func.sh 
cat func.log

echo "testing CaATPase case .. done"
cd ..
fi
#------------------------------------------------------------------------


echo
echo " ------ NORMA INSTALLATION TEST TERMINATED ------ "

echo "add the following line to your .bashrc or run the following command prior to using NORMA"
echo
echo "source "`pwd`"/software/setup.sh"
echo
