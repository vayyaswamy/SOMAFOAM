#!/bin/bash
#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.0
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     etc/settings.sh
#
# Description
#     Startup file for FOAM
#     Sourced from FOAM-??/etc/bashrc
#
#------------------------------------------------------------------------------

# prefix to PATH
_foamAddPath()
{
    while [ $# -ge 1 ]
    do
        export PATH=$1:$PATH
        shift
    done
}

# prefix to LD_LIBRARY_PATH
_foamAddLib()
{
    while [ $# -ge 1 ]
    do
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        if [ "$WM_ARCH_BASE" = "darwin" ]
        then
            # do NOT add the lib of MacPort as this might break programs
            if [ "$1" != "/opt/local/lib" ]
            then
                export DYLD_LIBRARY_PATH=$1:$DYLD_LIBRARY_PATH
            fi
        fi
        shift
    done
}

# Source files, possibly with some verbosity
# Yes, this is the same definition as in the file etc/bash
# We need that definition available for scripts sourcing
# settings.sh directly.
_foamSource()
{
   while [ $# -ge 1 ]
   do
      [ "$FOAM_VERBOSE" -a "$PS1" ] && echo "Sourcing: $1"
      . $1
      shift
   done
}

# location of the jobControl directory
export FOAM_JOB_DIR=$WM_PROJECT_INST_DIR

# wmake configuration
export WM_DIR=$WM_PROJECT_DIR/wmake
export WM_LINK_LANGUAGE=c++
export WM_OPTIONS=$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION
export PATH=$WM_DIR:$PATH

# base configuration
export FOAM_APP=$WM_PROJECT_DIR/applications
export FOAM_APPBIN=$WM_PROJECT_DIR/bin/$WM_OPTIONS
export FOAM_LIB=$WM_PROJECT_DIR/lib
export FOAM_LIBBIN=$WM_PROJECT_DIR/lib/$WM_OPTIONS
export FOAM_SRC=$WM_PROJECT_DIR/src

# shared site configuration - similar naming convention as ~FOAM expansion
export FOAM_SITE_DIR=$WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION
export FOAM_SITE_APPBIN=$FOAM_SITE_DIR/bin/$WM_OPTIONS
export FOAM_SITE_LIBBIN=$FOAM_SITE_DIR/lib/$WM_OPTIONS


# add FOAM scripts and wmake to the path
export PATH=$WM_DIR:$WM_PROJECT_DIR/bin:$PATH

_foamAddPath $FOAM_APPBIN $FOAM_SITE_APPBIN
_foamAddLib  $FOAM_LIBBIN $FOAM_SITE_LIBBIN


# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset compilerBin compilerLib

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compilerInstall = FOAM | System
#
# We can override the value of compilerInstall from prefs.sh
: ${compilerInstall:=System}

# Or we can force it right here
#compilerInstall=FOAM
compilerInstall=System

case "${compilerInstall}" in
FOAM)
    case "$WM_COMPILER" in
    Gcc)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.1.2/platforms/$WM_OPTIONS/etc/gmp-5.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.1.2/platforms/$WM_OPTIONS/etc/mpfr-3.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-1.0.1/platforms/$WM_OPTIONS/etc/mpc-1.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS/etc/gcc-4.6.4.sh
        ;;
    esac

    # Check that the compiler directory can be found
    if [ ! -d "$WM_COMPILER_DIR" ]
    then
        echo
        echo
    fi

    compilerBin=$WM_COMPILER_DIR/bin
    compilerLib=$WM_COMPILER_DIR/lib$WM_COMPILER_LIB_ARCH:$WM_COMPILER_DIR/lib
    ;;
esac

if [ -d "$compilerBin" ]
then
    _foamAddPath $compilerBin
    _foamAddLib  $compilerLib
fi

unset compilerBin compilerLib compilerInstall


if [ -z "$WM_CC" ]
then
    case "$WM_COMPILER" in
    Gcc*)
        export WM_CC='gcc'
        export WM_CXX='g++'
        ;;
    Icc)
        export WM_CC='icc'
        export WM_CXX='icpc'
        ;;
    esac
fi

# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH
mpi_version=unknown
case "$WM_MPLIB" in
OPENMPI)
    if [ ! -z $WM_THIRD_PARTY_USE_OPENMPI_188 ] && [ -e $WM_THIRD_PARTY_DIR/packages/openmpi-1.8.8/platforms/$WM_OPTIONS ]
        then
        mpi_version=openmpi-1.8.8
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using openmpi-1.8.8 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh
    fi

    # On Windows set mpi_version from value defined in bashrc.mingw:
    if [ "$WM_ARCH_BASE" == "mingw" ]
    then
        mpi_version=$MPI_VERSION_MINGW
    fi

    unset mpi_version
    ;;

SYSTEMOPENMPI)
    mpi_version=openmpi-system

    # make sure not the "old" mpi is used
    # Not sure if this is necessary anymore.
    # export OPAL_PREFIX=

    # Make sure OPENMPI_BIN_DIR is set and valid
    if [ -n "${OPENMPI_BIN_DIR}" ] && [ -d "${OPENMPI_BIN_DIR}" ]
    then
    # User defined value specified for OPENMPI_BIN_DIR
    #
    # WARNING:
    #          We assume this path specified by $OPENMPI_BIN_DIR is valid
    #          We assume the command mpicc is located somewhere under this path
    #          We assume the file mpi.h is located somewhere under this path
    #
    #          Otherwise, please double check your openmpi installation, you are
    #          probably missing the openmpi runtime and/or development packages
    #          available for your system.
    #
    _foamAddPath $OPENMPI_BIN_DIR
    else
    # Here, we assume your environment is already set for running
    # and developping with openmpi.
    #
    # Initialize OPENMPI_BIN_DIR using the path to mpicc
    export OPENMPI_BIN_DIR=$(dirname `which mpicc`)
    fi

    # Make sure OPENMPI_LIB_DIR is set
    if [ ! -n "${OPENMPI_LIB_DIR}" ]
    then
        # Initialize OPENMPI_LIB_DIR using the path to mpicc
        export OPENMPI_LIB_DIR="`mpicc --showme:libdirs`"
    fi

    # Make sure the dynamic libraries are accessible
    [   -n "${OPENMPI_LIB_DIR}" ]       && _foamAddLib $OPENMPI_LIB_DIR

    export MPI_HOME=`dirname $OPENMPI_BIN_DIR`
    export MPI_ARCH_PATH=$MPI_HOME
    export OPAL_PREFIX=$MPI_ARCH_PATH

    # We initialize the rest of the environment using mpicc --showme:
    [ ! -n "${OPENMPI_INCLUDE_DIR}" ]   && export OPENMPI_INCLUDE_DIR="`mpicc --showme:incdirs`"
    [ ! -n "${OPENMPI_COMPILE_FLAGS}" ] && export OPENMPI_COMPILE_FLAGS="`mpicc --showme:compile`"
    [ ! -n "${OPENMPI_LINK_FLAGS}" ]    && export OPENMPI_LINK_FLAGS="`mpicc --showme:link`"

    #
    # WARNING: We assume the file mpi.h will be available under the directories identified
    #          by the variable $OPENMPI_INCLUDE_DIR. Otherwise, please double check your
    #          system openmpi installation.

    # Set compilation flags here instead of in wmake/rules/../mplibSYSTEMOPENMPI
    export PINC="${OPENMPI_COMPILE_FLAGS}"
    export PLIBS="${OPENMPI_LINK_FLAGS}"

    # No longer needed, but we keep this as a reference, just in case...
    #libDir=`echo "$PLIBS" | sed -e 's/.*-L\([^ ]*\).*/\1/'`
    #_foamAddLib $libDir

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "  Environment variables defined for OpenMPI:"
        echo "    OPENMPI_BIN_DIR       : $OPENMPI_BIN_DIR"
        echo "    OPENMPI_LIB_DIR       : $OPENMPI_LIB_DIR"
        echo "    OPENMPI_INCLUDE_DIR   : $OPENMPI_INCLUDE_DIR"
        echo "    OPENMPI_COMPILE_FLAGS : $OPENMPI_COMPILE_FLAGS"
        echo "    OPENMPI_LINK_FLAGS    : $OPENMPI_LINK_FLAGS"
        echo ""
        echo "    MPI_HOME              : $MPI_HOME"
        echo "    MPI_ARCH_PATH         : $MPI_ARCH_PATH"
        echo "    OPAL_PREFIX           : $OPAL_PREFIX"
        echo "    PINC                  : $PINC"
        echo "    PLIBS                 : $PLIBS"
    fi

    unset mpi_version
    ;;

MPICH)
    mpi_version=mpich-1.2.4
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$mpi_version
    export MPI_ARCH_PATH=$MPI_HOME/platforms/$WM_OPTIONS
    export MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib

    unset mpi_version
    ;;

INTELMPI)
    # no trailing slash
    [ "${MPI_ROOT%/}" = "${MPI_ROOT}" ] || MPI_ROOT="${MPI_ROOT%/}"

    export FOAM_MPI="${MPI_ROOT##*/}"
    export MPI_ARCH_PATH=$MPI_ROOT

    if [ ! -d "$MPI_ROOT" -o -z "$MPI_ARCH_PATH" ]
    then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    MPI_ROOT not a valid mpt installation directory or ending" \
             " in a '/'." 1>&2
        echo "    Please set MPI_ROOT to the mpt installation directory." 1>&2
        echo "    MPI_ROOT currently set to '$MPI_ROOT'" 1>&2
    fi

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using INTEL MPI:" 1>&2
        echo "    MPI_ROOT : $MPI_ROOT" 1>&2
        echo "    FOAM_MPI : $FOAM_MPI" 1>&2
    fi

    _foamAddPath    $MPI_ARCH_PATH/intel64/bin
    _foamAddLib     $MPI_ARCH_PATH/intel64/lib

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/intelmpi
    ;;

*)
    ;;
esac

# Removed $FOAM_MPI_LIBBIN.  HJ, 8/Aug/2015


# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
minBufferSize=20000000

if [ "${MPI_BUFFER_SIZE:=$minBufferSize}" -lt $minBufferSize ]
then
    MPI_BUFFER_SIZE=$minBufferSize
fi
export MPI_BUFFER_SIZE

# Load Metis library
# ~~~~~~~~~~~~~~~~~~
[ -z "$METIS_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_METIS_510 ] && [ -e $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS/etc/metis-5.1.0.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    METIS_DIR is initialized to: $METIS_DIR"


# Load cmake
# ~~~~~~~~~~
[ -z "$CMAKE_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_CMAKE_322 ] && [ -e $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS/etc/cmake-3.2.2.sh
}

[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    CMAKE_DIR is initialized to: $CMAKE_DIR"

# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unset _foamAddPath _foamAddLib minBufferSize

# -----------------------------------------------------------------------------
