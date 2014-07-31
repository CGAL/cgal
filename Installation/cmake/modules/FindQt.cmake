#This file allows an easy add of Qt library, depending of Qt version used into CGAL compilation.
#New for Qt5 version !

if(${CGAL_Qt_version} )
	message("It seem that CGAL is not properly configured with Qt.")
else(${CGAL_Qt_version})

	find_package(${CGAL_Qt_version})
	
	if(QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
	endif(QT4_FOUND)
	
	#To the functions differences between Qt4 and Qt5.	
	include(QtChoice)
	
endif(${CGAL_Qt_version})