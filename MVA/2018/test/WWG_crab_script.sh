#this is not mean to be run locally
#
echo Check if TTY
if [ "`tty`" != "not a tty" ]; then
  echo "YOU SHOULD NOT RUN THIS IN INTERACTIVE, IT DELETES YOUR LOCAL FILES"
else

echo "ENV..................................."
env 
echo "VOMS"
voms-proxy-info -all
echo "CMSSW BASE, python path, pwd"
echo $CMSSW_BASE 
echo $PYTHON_PATH
echo $PWD 
rm -rf $CMSSW_BASE/lib/
rm -rf $CMSSW_BASE/src/
rm -rf $CMSSW_BASE/module/
rm -rf $CMSSW_BASE/python/
mv lib $CMSSW_BASE/lib
mv src $CMSSW_BASE/src
mv module $CMSSW_BASE/module
mv python $CMSSW_BASE/python
#mv ../WWG/* .

kind=""
mode=""
year=""
which_data=""
for i in "$@"
do
  case $i in
      kind=*)
      kind="${i#*=}"
      ;;
  esac
  case $i in
      mode=*)
      mode="${i#*=}"
      ;;
  esac
  case $i in
      year=*)
      year="${i#*=}"
      ;;
  esac
    case $i in
      which_data=*)
      which_data="${i#*=}"
      ;;
  esac
done

echo Found Proxy in: $X509_USER_PROXY
#if [ $isdata == "1" ]; then
python WWG_postproc.py -k $kind -m $mode -y $year -w $which_data
#else
#  if [ $iswwa == "1" ]; then
#      python whjj_postproc.py -k $kind -m -y $year
#  else
#     python whjj_postproc.py -y $year
#  fi

fi
