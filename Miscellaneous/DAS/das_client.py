# !/bin/sh
# VERSION:slc7_aarch64_gcc700/v03.01.00-omkpbe3
# Clean-up CMSSW environment
eval `scram unsetenv -sh 2>/dev/null`
# Sourcing dasclient environment
SHARED_ARCH=`cmsos`
LATEST_VERSION=`cd /cvmfs/cms.cern.ch; ls ${SHARED_ARCH}_*/cms/das_client/v*/etc/profile.d/init.sh | sed 's|.*/cms/das_client/||' | sort | tail -1`
DAS_ENV=`ls /cvmfs/cms.cern.ch/${SHARED_ARCH}_*/cms/das_client/${LATEST_VERSION} | sort | tail -1`
source $DAS_ENV
if [ $# == 0 ] || [ "$1" == "--help" ] || [ "$1" == "-help" ]
then
    $DAS_CLIENT_ROOT/bin/das_client.py --help | sed 's/das_client.py/das_client/'
else
    $DAS_CLIENT_ROOT/bin/das_client.py "$@"
fi
