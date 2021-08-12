#
# Example environment setup script (called by example_config.sh)
#
#!/bin/bash

export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH


export ROOTSYS=/vols/t2k/users/pjd12/analysiswork/t2kreweightthings/psychedir/ROOT/v5r34p34n00/Linux-x86_64


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

# For Mac's
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
