# Load metis-5.1.0 libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export METIS_DIR=$WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS
export METIS_BIN_DIR=$METIS_DIR/bin
export METIS_LIB_DIR=$METIS_DIR/lib
export METIS_INCLUDE_DIR=$METIS_DIR/include

# Enable access to the package runtime applications and libraries
[ -d $METIS_BIN_DIR ] && _foamAddPath $METIS_BIN_DIR
[ -d $METIS_LIB_DIR ] && _foamAddLib  $METIS_LIB_DIR

