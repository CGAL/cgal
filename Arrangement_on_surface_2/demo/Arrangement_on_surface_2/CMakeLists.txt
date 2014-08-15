# This is the CMake script for compiling a CGAL application.

project( Arrangement_2 ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set( QT_USE_QTOPENGL  TRUE )

find_package(CGAL REQUIRED COMPONENTS Core Qt4or5)

include(${CGAL_USE_FILE})
	
find_package(Qt4or5)

if ( CGAL_FOUND AND CGAL_Core_FOUND AND CGAL_${QT_VERSION}_FOUND AND QT${QT_VERSION_USED}_FOUND )

  include_directories( ./ )
  # Arrangement package includes
  include_directories( BEFORE ../../include )

 
  qt_wrap_ui( arrangement_2_uis
    ArrangementDemoWindow.ui
    NewTabDialog.ui
    OverlayDialog.ui
    ArrangementDemoPropertiesDialog.ui
  )
  
  if(QT4)
  qt_wrap_cpp( arrangement_2_mocs
    ArrangementDemoWindow.h
    ArrangementDemoTab.h
    Callback.h
#    NewTabDialog.h
    OverlayDialog.h
    ArrangementDemoPropertiesDialog.h
    ColorItemEditor.h
    DeleteCurveModeItemEditor.h
    PropertyValueDelegate.h
#    PropertyValueDelegate.cpp
  )
  endif(QT4)

  qt_add_resources( arrangement_2_resources
    ArrangementDemoWindow.qrc
  )
  add_executable( arrangement_2
    arrangement_2.cpp
    ArrangementGraphicsItem.cpp
    ArrangementDemoWindow.cpp
    ArrangementDemoTab.cpp
    ArrangementDemoGraphicsView.cpp
    Callback.cpp
    VerticalRayShootCallback.cpp
    EnvelopeCallback.cpp
    SplitEdgeCallback.cpp
    FillFaceCallback.cpp
    GraphicsViewSegmentInput.cpp
    GraphicsViewCurveInput.cpp
    Utils.cpp
    NewTabDialog.cpp
    OverlayDialog.cpp
    ArrangementDemoPropertiesDialog.cpp
    ColorItemEditor.cpp
    PropertyValueDelegate.cpp
    DeleteCurveMode.cpp
    DeleteCurveModeItemEditor.cpp
    PointsGraphicsItem.cpp
    VerticalRayGraphicsItem.cpp
   ${arrangement_2_mocs}
   ${arrangement_2_uis}
   ${arrangement_2_resources}
  )

  target_link_libraries( arrangement_2
    ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
    ${QT_LIBRARIES}
  )

else( CGAL_FOUND AND CGAL_Core_FOUND AND CGAL_${QT_VERSION}_FOUND AND QT${QT_VERSION_USED}_FOUND )
    
  set(ARRANGEMENT_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(ARRANGEMENT_MISSING_DEPS "the CGAL library, ${ARRANGEMENT_MISSING_DEPS}")
  endif()
    
  if(NOT CGAL_Core_FOUND)
    set(ARRANGEMENT_MISSING_DEPS "the CGAL Core library, ${ARRANGEMENT_MISSING_DEPS}")
  endif()
  
  if(NOT CGAL_${QT_VERSION}_FOUND)
    set(ARRANGEMENT_MISSING_DEPS "the CGAL Qt${QT_VERSION_USED} library, ${ARRANGEMENT_MISSING_DEPS}")
  endif()
  
  if(NOT QT${QT_VERSION_USED}_FOUND)
    set(ARRANGEMENT_MISSING_DEPS "Qt${QT_VERSION_USED}, ${ARRANGEMENT_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: The Arrangemnt on surface_2 demonstration requires ${ARRANGEMENT_MISSING_DEPS}and will not be compiled.")

endif()

# add_subdirectory( tests )
