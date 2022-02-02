# Load metis-5.1.0 libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv METIS_DIR $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS
setenv METIS_BIN_DIR $METIS_DIR/bin
setenv METIS_LIB_DIR $METIS_DIR/lib
setenv METIS_INCLUDE_DIR $METIS_DIR/include

if ( -e $METIS_BIN_DIR ) then
    _foamAddPath $METIS_BIN_DIR
endif

if ( -e $METIS_LIB_DIR ) then
    _foamAddLib $METIS_LIB_DIR
endif
