# -----------------------------------------------------------------------------
#  REST_SMFI | Aetup_REST_SMFI.sh
#
#  Setup enviorment script
#   * Author: Austin McDonald
#   * Creation date: Oct 2019
# -----------------------------------------------------------------------------


printf '************************************* \n'
message="Defining python paths"; for ((i=0; i<${#message}; i++)); do echo "after 75" | tclsh; printf "${message:$i:1}"; done; echo;
sleep 1

export SMFIDIR=$(pwd)

export COREDIR=$SMFIDIR/Core
export CYTHONDIR=$SMFIDIR/Cython

export PYTHONPATH=$COREDIR:$PYTHONPATH
export PYTHONPATH=$CYTHONDIR:$PYTHONPATH
export PATH=$PATH:$CYTHONDIR

#cd $CYTHONDIR
#bash buils.sh
#cd $SMFIDIR

printf '************************************* \n'

message="Paths defined!"; for ((i=0; i<${#message}; i++)); do echo "after 75" | tclsh; printf "${message:$i:1}"; done; echo;

printf '************************************* \n'
