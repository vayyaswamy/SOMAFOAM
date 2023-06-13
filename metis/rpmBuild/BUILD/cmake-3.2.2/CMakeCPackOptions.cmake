# This file is configured at cmake time, and loaded at cpack time.
# To pass variables to cpack from cmake, they must be configured
# in this file.

if(CPACK_GENERATOR MATCHES "NSIS")
  set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")

  # set the install/unistall icon used for the installer itself
  # There is a bug in NSI that does not handle full unix paths properly.
  set(CPACK_NSIS_MUI_ICON "/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Utilities/Release\\CMakeLogo.ico")
  set(CPACK_NSIS_MUI_UNIICON "/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Utilities/Release\\CMakeLogo.ico")
  # set the package header icon for MUI
  set(CPACK_PACKAGE_ICON "/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Utilities/Release\\CMakeInstall.bmp")
  # tell cpack to create links to the doc files
  set(CPACK_NSIS_MENU_LINKS
    "doc/cmake-3.2/html/index.html" "CMake Documentation"
    "http://www.cmake.org" "CMake Web Site"
    )
  # Use the icon from cmake-gui for add-remove programs
  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\cmake-gui.exe")

  set(CPACK_NSIS_PACKAGE_NAME "CMake 3.2.2")
  set(CPACK_NSIS_DISPLAY_NAME "CMake 3.2.2, a cross-platform, open-source build system")
  set(CPACK_NSIS_HELP_LINK "http://www.cmake.org")
  set(CPACK_NSIS_URL_INFO_ABOUT "http://www.kitware.com")
  set(CPACK_NSIS_CONTACT cmake@cmake.org)
  set(CPACK_NSIS_MODIFY_PATH ON)
endif()

# include the cpack options for qt dialog if they exist
# they might not if qt was not enabled for the build
include("/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Source/QtDialog/QtDialogCPack.cmake" OPTIONAL)

if(CPACK_GENERATOR MATCHES "IFW")
  # Installer configuration
  set(CPACK_IFW_PACKAGE_TITLE "CMake Build Tool")
  set(CPACK_IFW_PRODUCT_URL "http://www.cmake.org")
  
  set(CPACK_IFW_PACKAGE_WINDOW_ICON
    "/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Source/QtDialog/CMakeSetup128.png")
  # Package configuration group
  set(CPACK_IFW_PACKAGE_GROUP CMake)
  # Group configuration
  set(CPACK_COMPONENT_GROUP_CMAKE_DISPLAY_NAME
    "CMake")
  set(CPACK_COMPONENT_GROUP_CMAKE_DESCRIPTION
    "CMake is a build tool")
  # IFW group configuration
  set(CPACK_IFW_COMPONENT_GROUP_CMAKE_VERSION
    "3.2.2")
  set(CPACK_IFW_COMPONENT_GROUP_CMAKE_LICENSES
    "CMake Copyright" "/data/vayyaswamy/ltps/metis/rpmBuild/BUILD/cmake-3.2.2/Copyright.txt")
  
endif()

if(CPACK_GENERATOR MATCHES "CygwinSource")
  # when packaging source make sure the .build directory is not included
    set(CPACK_SOURCE_IGNORE_FILES
      "/CVS/" "/\\.build/" "/\\.svn/" "\\.swp$" "\\.#" "/#" "~$")
endif()

if("${CPACK_GENERATOR}" STREQUAL "PackageMaker")
  if(CMAKE_PACKAGE_QTGUI)
    set(CPACK_PACKAGE_DEFAULT_LOCATION "/Applications")
  else()
    set(CPACK_PACKAGE_DEFAULT_LOCATION "/usr")
  endif()
endif()

if("${CPACK_GENERATOR}" STREQUAL "WIX")
  # Reset CPACK_PACKAGE_VERSION to deal with WiX restriction.
  # But the file names still use the full CMake_VERSION value:
  set(CPACK_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}-3.2.2-${CPACK_SYSTEM_NAME}")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}-3.2.2-Source")

  if(NOT CPACK_WIX_SIZEOF_VOID_P)
    set(CPACK_WIX_SIZEOF_VOID_P "8")
  endif()

  set(CPACK_PACKAGE_VERSION
    "3.2")
  # WIX installers require at most a 4 component version number, where
  # each component is an integer between 0 and 65534 inclusive
  set(patch "2")
  if(patch MATCHES "^[0-9]+$" AND patch LESS 65535)
    set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}.${patch}")
  endif()
endif()
