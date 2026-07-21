# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin
EDITOR=emacs
export EDITOR
export PATH

# ENV variables that need to be defined for PAWS make to work:
export MATLABROOT=/gpfs/apps/x86_64-rhel6/matlab/R2012b
# These should normally stay the same:
export LD_LIBRARY_PATH=$MATLABROOT/bin/glnxa64:$MATLABROOT/sys/os/glnxa64:/lib64:$LD_LIBRARY_PATH
export MLIBS=$MATLABROOT/bin/glnxa64/
export INCLUDE_PATH=$MATLABROOT/extern/include
ulimit -s unlimited