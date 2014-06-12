message("Searching Qt5 modules.")

if(WIN32)
	message("Qt5 on Windows needs Windows SDK.")
	
	set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:\\Program Files (x86)\\Windows Kits\\8.1\\Lib\\winv6.3\\um\\x64")
	
endif()

# Qt modules
FOREACH(module Core GUI OpenGL Multimedia  
		Network QML Quick SQL Test WebKit 
		Widgets D-Bus Graphical_Effects ImageFormats 
		MacExtras NFC Positioning PrintSupport Declarative 
		Script Sensors SerialPort SVG WebSockets WindowsExtras
		X11Extras Xml XMLPatterns Designer Help UITools)

  STRING(TOUPPER ${module} component)

  IF (QT_USE_QT${component})
    	FIND_PACKAGE(Qt5${module} REQUIRED)
    	IF (Qt5${module}_FOUND)
      	message(STATUS "Qt5${module} found.")
      	SET(QT_INCLUDE_DIR ${QT_INCLUDE_DIR} ${Qt5${module}_INCLUDE_DIRS})
      	SET(QT_LIBRARIES ${QT_LIBRARIES} ${Qt5${module}_LIBRARIES})
      	SET(QT_DEFINITIONS ${QT_DEFINITIONS} ${Qt5${module}_DEFINITIONS})
    	ELSE (Qt5${module}_FOUND)
      		MESSAGE("Qt5 ${module} library not found.")
    	ENDIF (Qt5${module}_FOUND)
  ENDIF (QT_USE_QT${component})
  
ENDFOREACH(module)

SET(QT5 TRUE)
set(CMAKE_AUTOMOC ON)
