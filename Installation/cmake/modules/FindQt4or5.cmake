#This file allows an easy add of Qt library with CGAL.
#New for Qt5 version !
#
#
#find_package(Qt4or5) works with the variable ${CGAL_QT_VERSION}, containing the value of Qt version using.
#
#So, the variable ${CGAL_Qt_version} must be initialised before use find_package(Qt4or5) if auto configuration wanted. 
#
#At the end of the process, QT_VERSION_USED is set to the Qt proper version number.
#

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

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
