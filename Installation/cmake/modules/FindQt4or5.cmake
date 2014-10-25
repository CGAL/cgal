#
#=================================================FindQt4or5.cmake=================================================
#
#	FindQt4or5.cmake is a new way to find package of Qt, with the possibility of make a choice between Qt4 or Qt5.
#
#
#	Several variables are linked to this file :
#  
#  - QT_VERSION_USED
#  - QT_VERSION
#  - USE_QT_VERSION
#  - CGAL_QT_VERSION
#
#
#==================================================QT_VERSION_USED==================================================

#	The variable QT_VERSION_USED is set when an Qt version is asked into FindQt4or5.cmake even if this Qt version 
# is not found.
#
#   Hence, if Qt4 asked,  QT_VERSION_USED value is set to 4 and to 5 if Qt5 asked.
#
#	Moreover, when an Qt version ( 4 or 5 ) is found, the followings flags are set :
#   
#  - QT4_FOUND => Qt4 has been loaded
#  - QT5_FOUND => Qt5 has been loaded
#
#	So, into CMakeList, the typical use of these flag could be : 
#
#   IF ( QT${QT_VERSION_USED)_FOUND ) //Works for QT4_FOUND and QT5_FOUND
#   ...
#
#
#=====================================================QT_VERSION====================================================
#
#	The value of QT_VERSION variable depends on the QT_VERSION_USED value. Indeed, if QT_VERSION_USED value is 
# 5, so the QT_VERSION will be "Qt5" and Qt4 if QT_VERSION_USED is 4.
#
#
#
#==================================================USE_QT_VERSION==================================================
#
# When find_package(Qt4or5) is used, two possibility are offered
#
#
# 1) If the variable USE_QT_VERSION isn't set, an drop-down menu for Qt version choice appear on CMake-gui
#
# 2) If the variable USE_QT_VERSION is set, so the user force the Qt version ( 4 or 5 ) to use and no choice 
#    is proposed into CMake-gui
#
#	Of course, during the process, QT_VERSION_USED is set to the Qt proper version number.
#
#
#==================================================CGAL_QT_VERSION==================================================
#
#	The last variable used into FindQt4or5.cmake is CGAL_QT_VERSION, an variable set during the process 
# of USE_CGAL.cmake.
#
#	Actually, when find_package(CGAL COMPONENT Qt4or5) is called, the variable ${CGAL_QT_VERSION} is set and when 
# find_package(Qt4or5) is called (after), the Qt version is set to the proper Qt version 
# chosen for CGAL.
#
# 	For instance, if libCGAL_Qt5 chosen during find_package(CGAL COMPONENTS Qt4or5) process, Qt5 will be automatically 
# loaded when find_package(Qt4or5) calling.
#


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

#Get the previous state of QT_VERSION.
cache_get(QT_VERSION_mem)
set(OLD_QT_VERSION ${QT_VERSION_mem})


#Set of a temporary variable, because the if existence checking seem to not work without it.
set (USE_QT_VERSION_temp ${USE_QT_VERSION})

if(USE_QT_VERSION_temp)
  set(QT_VERSION_CHOICE FALSE)
  set(QT_VERSION Qt${USE_QT_VERSION})
else()
  set(QT_VERSION_CHOICE TRUE)
endif()

set (OLDCGALQt_version_temp ${OLD_CGAL_QT_VERSION})
set (CGALQt_version_temp ${CGAL_QT_VERSION})

if (OLDCGALQt_version_temp AND CGALQt_version_temp)
  if(NOT (${CGAL_QT_VERSION} STREQUAL ${OLD_CGAL_QT_VERSION}) )
    UNSET(QT_CHOICE CACHE)
  endif()
endif()

#Same as before. 
set (CGAL_QT_VERSION_temp ${CGAL_QT_VERSION})

if(NOT CGAL_QT_VERSION_temp)
  set(CGAL_QT_VERSION "5")
else()
  message(STATUS "Qt configuration : libCGAL_Qt${CGAL_QT_VERSION} has been asked.")
endif()

if(${QT_VERSION_CHOICE})
  cache_set(OLD_CGAL_QT_VERSION ${CGAL_QT_VERSION})

  set (QT_CHOICE ${CGAL_QT_VERSION} CACHE STRING "Choice of Qt version for find_package(Qt4or5).")
  SET_PROPERTY(CACHE QT_CHOICE PROPERTY STRINGS 4 5)
  set(QT_VERSION Qt${QT_CHOICE})

  #Save the current version of Qt considered.
  cache_set(QT_VERSION_mem ${QT_VERSION})
endif()

#Check if we switch Qt version (from 4 to 5 or 5 to 4)
if( OLD_QT_VERSION )
  if(NOT ${OLD_QT_VERSION} STREQUAL ${QT_VERSION})
    message("Switch from ${OLD_QT_VERSION} to ${QT_VERSION}")
    message("Think to change the version of externals libraries that depending on Qt and to clean the build directory before make the project.")
    set(QT_QMAKE_CHANGED TRUE)
  endif()
endif()

if(${QT_VERSION} STREQUAL "Qt4")
  UNSET(QT4} CACHE)
  UNSET(QT4_FOUND CACHE)
  UNSET(QT_VERSION_USED CACHE)

  #We say that we want the version 4 of the Qt library.
  set(QT_VERSION_USED 4)
endif()

find_package(${QT_VERSION})

if(${QT_VERSION} STREQUAL "Qt4" AND QT4_FOUND)
  include(${QT_USE_FILE})
  set(QT4 TRUE)

  UNSET(QT4 CACHE)
  UNSET(QT5_fOUND CACHE)

  message("Qt4 found")
endif()

#To the functions differences between Qt4 and Qt5.
include(QtChoice)
