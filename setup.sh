#
# Example environment setup script (called by example_config.sh)
#
#!/bin/bash


#SETUP CMAKE AND ROO FROM CVMFS
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

spack load python@3.9.15
echo “root”
spack load root@6.28.12
echo “cmake”
spack load cmake@3.27.7
echo “gcc”
spack load gcc@12.2.0%gcc@11.4.1
spack load fife-utils@3.7.0
spack load metacat@4.0.0
spack load rucio-clients@33.3.0
spack load sam-web-client@3.4%gcc@12.2.0
spack load r-m-dd-config@1.0 experiment=dune
export IFDH_CP_MAXRETRIES=0\0\0\0\0

export CXX=`which g++` # this might be specific for Fermilab?
export CC=`which gcc` # this might be specific for Fermilab?

export LD_LIBRARY_PATH=$PWD/src_lib:$LD_LIBRARY_PATH


# GENIE and dependencies
# Only needed with ./configure --enable-genie
#export GENIE=/path/to/genie
#export PYTHIA6_LIB=/path/to/pythia6/lib
#export LIBXML_LIB=/path/to/libxml/lib
#export LOG4CPP_LIB=/path/to/log4cpp/lib
#export LHAPDF_LIB=/path/to/lhapdf/lib
#export LD_LIBRARY_PATH=$LOG4CPP_LIB:$LIBXML_LIB:$LHAPDF_LIB:$PYTHIA6_LIB:$ROOTSYS/lib:$GENIE/lib:$LD_LIBRARY_PATH;
#export LIBXML_INC=/path/to/libxml2/include
#export LOG4CPP_INC=/path/to/log4cpp/include
#export LHAPDF_INC=/path/to/lhapdf/include/
#export LHAPATH=/path/to/lhapdf/PDFsets
#export PATH=$GENIE/bin:$PATH;
#
# # For Mac's
# export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
