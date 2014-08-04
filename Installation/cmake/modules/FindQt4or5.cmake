#This file allows an easy add of Qt library with CGAL.
#New for Qt5 version !
#
#
#find_package(Qt4or5) works with the variable ${CGAL_QT_VERSION}, containing the value of Qt version using.
#
#So, the variable ${CGAL_Qt_version} must be initialised before use find_package(Qt4or5) if auto configuration wanted. 
#
#At the end of the process, USE_QT_VERSION is set to the Qt proper version number.
#


#Set of a temporary variable, because the if existence checking seem to not work without it.
set (Qt_version_temp ${Qt_version})

if(Qt_version_temp)
	set(QT_VERSION_CHOICE FALSE)
else()
	set(QT_VERSION_CHOICE TRUE)
endif()

set (OLDCGALQt_version_temp ${OLD_CGAL_QT_VERSION})
set (CGALQt_version_temp ${CGAL_QT_VERSION})

if (OLDCGALQt_version_temp AND CGALQt_version_temp)
	if(NOT (${CGAL_QT_VERSION} STREQUAL ${OLD_CGAL_QT_VERSION}) )
		UNSET(QT_VERSION CACHE)
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
	
		set (QT_VERSION ${CGAL_QT_VERSION} CACHE STRING "Choice of Qt version for find_package(Qt4or5).")
		SET_PROPERTY(CACHE QT_VERSION PROPERTY STRINGS 4 5)
		set(Qt_version Qt${QT_VERSION})
	endif()
	
	if(${Qt_version} STREQUAL "Qt4")
		UNSET(QT4} CACHE)
		UNSET(QT4_FOUND CACHE)
		UNSET(USE_QT_VERSION CACHE)
		
		#We say that we want the version 4 of the Qt library.
		set(USE_QT_VERSION 4)
    endif()
	
	find_package(${Qt_version})
	
	if(${Qt_version} STREQUAL "Qt4" AND QT4_FOUND)
		include(${QT_USE_FILE})
		set(QT4 TRUE)

		UNSET(QT4 CACHE)
		UNSET(QT5_fOUND CACHE)

		message("Qt4 found")
	endif()
	
	#To the functions differences between Qt4 and Qt5.	
	include(QtChoice)
