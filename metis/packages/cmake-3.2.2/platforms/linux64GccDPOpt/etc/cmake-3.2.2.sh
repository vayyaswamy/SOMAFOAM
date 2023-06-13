# Load cmake-3.2.2 libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export CMAKE_DIR=$WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS
export CMAKE_BIN_DIR=$CMAKE_DIR/bin

# Enable access to the runtime package applications
[ -d $CMAKE_BIN_DIR ] && _foamAddPath $CMAKE_BIN_DIR
