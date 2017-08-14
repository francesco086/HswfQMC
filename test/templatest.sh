#!/bin/bash
# First argument is expected to be the HswfQMC executable's full path 

TEST_PATH=$(pwd)

if [ "$1" == "" ]
  then
    HQMC_PATH=${TEST_PATH}/../HswfQMC_exe
  else
    HQMC_PATH=$1
fi

RUN_PATH=${TEST_PATH}/run/templatest.run
CMP_PATH=${TEST_PATH}/cmp/templatest.cmp
INP_PATH=${TEST_PATH}/inp/templatest.inp

rm -rv $RUN_PATH
cp -r $INP_PATH $RUN_PATH

cd $RUN_PATH
$HQMC_PATH > templatest.out

cd ../
diff -r $CMP_PATH $RUN_PATH > templatest.diff
if [ "$(wc -c templatest.diff)" == "0 templatest.diff" ]
  then
    echo "Test finished succesfully!"
  else
    echo "Test produces unexpected output:"
    cat templatest.diff
fi

cd $TEST_PATH
