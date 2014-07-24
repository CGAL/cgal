message("Searching Qt5 modules.")

if(WIN32)
	message("Qt5 on Windows needs Windows SDK.")
	
	set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:\\Program Files (x86)\\Windows Kits\\8.1\\Lib\\winv6.3\\um\\x64")
	
endif()

if (QT_USE_QTMAIN OR NOT component)
	find_package(Qt5Core)
endif(QT_USE_QTMAIN OR NOT component)

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

  #######################################
  #
  #       Check the executables of Qt 
  #          ( moc, uic, rcc )
  #         Same as Qt4 version  
  #
  #######################################


  IF(QT_QMAKE_CHANGED)
    SET(QT_UIC_EXECUTABLE NOTFOUND)
    SET(QT_MOC_EXECUTABLE NOTFOUND)
    SET(QT_UIC3_EXECUTABLE NOTFOUND)
    SET(QT_RCC_EXECUTABLE NOTFOUND)
    SET(QT_DBUSCPP2XML_EXECUTABLE NOTFOUND)
    SET(QT_DBUSXML2CPP_EXECUTABLE NOTFOUND)
    SET(QT_LUPDATE_EXECUTABLE NOTFOUND)
    SET(QT_LRELEASE_EXECUTABLE NOTFOUND)
    SET(QT_QCOLLECTIONGENERATOR_EXECUTABLE NOTFOUND)
    SET(QT_DESIGNER_EXECUTABLE NOTFOUND)
    SET(QT_LINGUIST_EXECUTABLE NOTFOUND)
  ENDIF(QT_QMAKE_CHANGED)
  
  FIND_PROGRAM(QT_MOC_EXECUTABLE
    NAMES moc-qt5 moc
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_UIC_EXECUTABLE
    NAMES uic-qt5 uic
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_UIC3_EXECUTABLE
    NAMES uic3
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_RCC_EXECUTABLE 
    NAMES rcc
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_DBUSCPP2XML_EXECUTABLE 
    NAMES qdbuscpp2xml
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_DBUSXML2CPP_EXECUTABLE 
    NAMES qdbusxml2cpp
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_LUPDATE_EXECUTABLE
    NAMES lupdate-qt5 lupdate
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_LRELEASE_EXECUTABLE
    NAMES lrelease-qt5 lrelease
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_QCOLLECTIONGENERATOR_EXECUTABLE
    NAMES qcollectiongenerator-qt5 qcollectiongenerator
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_DESIGNER_EXECUTABLE
    NAMES designer-qt5 designer
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  FIND_PROGRAM(QT_LINGUIST_EXECUTABLE
    NAMES linguist-qt5 linguist
    PATHS ${QT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  IF (QT_MOC_EXECUTABLE)
     SET(QT_WRAP_CPP "YES")
  ENDIF (QT_MOC_EXECUTABLE)

  IF (QT_UIC_EXECUTABLE)
     SET(QT_WRAP_UI "YES")
  ENDIF (QT_UIC_EXECUTABLE)



  MARK_AS_ADVANCED( QT_UIC_EXECUTABLE QT_UIC3_EXECUTABLE QT_MOC_EXECUTABLE
    QT_RCC_EXECUTABLE QT_DBUSXML2CPP_EXECUTABLE QT_DBUSCPP2XML_EXECUTABLE
    QT_LUPDATE_EXECUTABLE QT_LRELEASE_EXECUTABLE QT_QCOLLECTIONGENERATOR_EXECUTABLE
    QT_DESIGNER_EXECUTABLE QT_LINGUIST_EXECUTABLE)


set(QT5 TRUE)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("End of searching Qt5 modules.")
