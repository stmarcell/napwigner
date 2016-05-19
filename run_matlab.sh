#!/bin/bash


# Matlab setting
PREFDIR=/home/umat/.matlab/R2013a/users/stippinger/
# Matlab commands
PRECOMMAND="addpath(genpath('~/marcell/_DataHigh1.2/'),'~/marcell/napwigner/');"
ERRHANDLING="fprintf('ERR: %s\\n',ME.identifier); fprintf('%s\\n',ME.message); disp(struct2table(ME.stack)); celldisp(ME.cause);"

# Run simulation
MATLAB_PREFDIR=$PREFDIR matlab -nodisplay -r "$PRECOMMAND;"
#MATLAB_PREFDIR=$PREFDIR matlab -nodisplay -r "$PRECOMMAND; try; $mfunction; catch ME; $ERRHANDLING; end; exit;"

