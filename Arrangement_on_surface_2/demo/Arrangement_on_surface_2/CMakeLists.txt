# This is the CMake script for compiling a CGAL application.

project( Arrangement_2 ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL COMPONENTS Core Qt5 )
include( ${CGAL_USE_FILE} )

find_package( Qt5 COMPONENTS Script Widgets )

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)
  include_directories( ./ )
  # Arrangement package includes
  include_directories( BEFORE ../../include )

  qt5_wrap_ui( arrangement_2_uis
    ArrangementDemoWindow.ui
    NewTabDialog.ui
    OverlayDialog.ui
    ArrangementDemoPropertiesDialog.ui
  )
  qt5_wrap_cpp( arrangement_2_mocs
    ArrangementDemoWindow.h
    ArrangementDemoTab.h
    Callback.h
#    NewTabDialog.h
    OverlayDialog.h
    ArrangementDemoPropertiesDialog.h
    ColorItemEditor.h
    DeleteCurveModeItemEditor.h
    PropertyValueDelegate.h
    #PropertyValueDelegate.cpp
  )
  qt5_add_resources( arrangement_2_resources
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
  qt5_use_modules(arrangement_2 Widgets Script)
  target_link_libraries( arrangement_2
    ${CGAL_LIBRARIES}
    ${CGAL_3RD_PARTY_LIBRARIES}
    ${QT_LIBRARIES}
  )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS arrangement_2 )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()

# add_subdirectory( tests )
