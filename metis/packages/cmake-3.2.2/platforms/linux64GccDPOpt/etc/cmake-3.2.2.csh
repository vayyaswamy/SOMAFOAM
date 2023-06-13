# Load cmake-3.2.2 libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv CMAKE_DIR $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS
setenv CMAKE_BIN_DIR $CMAKE_DIR/bin

if ( -e $CMAKE_BIN_DIR ) then
    _foamAddPath $CMAKE_BIN_DIR
endif
