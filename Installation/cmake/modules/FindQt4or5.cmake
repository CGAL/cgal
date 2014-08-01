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


if (${OLD_CGAL_QT_VERSION})
	if(NOT (${CGAL_QT_VERSION} STREQUAL ${OLD_CGAL_QT_VERSION}) )
		UNSET(QT_VERSION CACHE)
	endif()
endif()
 

if(NOT ${CGAL_QT_VERSION})
	set(CGAL_QT_VERSION "5")
else(NOT ${CGAL_QT_VERSION})
    message(STATUS "Qt configuration : libCGAL_Qt${CGAL_QT_VERSION} has been loaded.")
endif(NOT ${CGAL_QT_VERSION})
    
	cache_set(OLD_CGAL_QT_VERSION ${CGAL_QT_VERSION})
	
	set (QT_VERSION ${CGAL_QT_VERSION} CACHE STRING "Choice of Qt version for find_package(Qt4or5).")
	SET_PROPERTY(CACHE QT_VERSION PROPERTY STRINGS 4 5)
	set(Qt_version Qt${QT_VERSION})
	
	find_package(${Qt_version})
	
	if(QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
		set(USE_QT_VERSION 4)
	endif(QT4_FOUND)
	
	#To the functions differences between Qt4 and Qt5.	
	include(QtChoice)