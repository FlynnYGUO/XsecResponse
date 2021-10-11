Create a file full of splines that are the response functions in each bin of true/reco energy for the Nue analysis.  

To Use (from top directory):
make
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH 
bin/make_xsec_response \
   -w /path/to/caffile.root \
   -m /path/to/wgtfile.root \
   -o path/to/desired/outputfile.root
   -selec {nue or numu} \

This will use the default templates in the input directory.  You can specify alternative templates, to see the arguments run the above executable with no command line arguments.  There is a small example script for making a nue reconstructed bin template:

inputs/make_template.C
