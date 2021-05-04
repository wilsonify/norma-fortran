#FILE:    install.sh
#AUTHOR:  Karsten Suhre
#DATE:    Tue Apr 18 16:55:10 CEST 2006
#PURPOSE: install NORMA
#BUGS:    
#MODIF:   

if [ ! -e ./software ] ; then

   echo "ERROR: no subdirectory ./software found"
  exit 1

fi

# generate setup files
echo "generating setup scripts .."
cat > ./software/setup.sh <<EOF1
# NORMA setup file for bash
# set the PATH to the software directory
export PATH=`pwd`/software:\$PATH

# set PATH and BIN variable for  URO
# if you are using your own URO version, adapt these lines
# to your local URO installation
export URO=`pwd`/software/URO
export BIN="bin"

# on some systems, the following variables need to be set to en_US
# in order to prevent real numbers to be treated using commata rather dots,
# otherwise this would lead to processing errors in some URO awk scripts
export LC_ALL="en_US"
export LANG="en_US"

EOF1

cat > ./software/setup.csh <<EOF2
# NORMA setup file for csh
# set the PATH to the software directory
setenv PATH `pwd`/software:\$PATH

# set PATH and BIN variable for  URO
# if you are using your own URO version, adapt these lines
# to your local URO installation
setenv URO `pwd`/software/URO
setenv BIN "bin"

# on some systems, the following variables need to be set to en_US
# in order to prevent real numbers to be treated using commata rather dots,
# otherwise this would lead to processing errors in some URO awk scripts
setenv LC_ALL "en_US"
setenv LANG "en_US"

EOF2
echo "generating setup script .. done"

# compile NORMA code (if required)
cd software
make
cd ..

# run tests (and complete installation)
./test.sh NORMA

