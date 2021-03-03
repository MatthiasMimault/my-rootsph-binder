name=DBG
casename=Dbg-arallelDiv_1
dirout=${name}_out
options=s
n_avg=4

# "executables" are renamed and called from their directory
dirbin=bin/windows
pyStats="${dirbin}/pystats_main.py"


# Library path must be indicated properly
current=$(pwd)
cd $dirbin
path_so=$(pwd)
cd $current
export LD_LIBRARY_PATH=$path_so
diroutdata=${dirout}/data;


# CODES are executed according the selected parameters of execution in this testcase
errcode=0

dircsv=${dirout}/csv;
dirimg=${dirout}/img;

if [ $errcode -eq 0 ]; then
  # Call stats to img script
  python $pyStats $casename $dircsv $dirimg $options $n_avg
  errcode=$?
fi


if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi

read -n1 -r -p "Press any key to continue..." key
echo
