#!/bin/bash

cd ${0%/*} || exit 1    # run from this directory

source "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/etc/bashrc

echo "*---------------------------*LTPS Installation Initiated*--------------------------------*"
echo

# wmake is required for subsequent targets
( cd wmake/src && make )

# build ThirdParty sources
( cd metis && ./AllMake )

. $WM_PROJECT_DIR/etc/settings.sh

# build all libraries and applications
src/Allwmake

applications/Allwmake

echo
echo "*--------------------------*LTPS Installation Finished*----------------------------------*"
