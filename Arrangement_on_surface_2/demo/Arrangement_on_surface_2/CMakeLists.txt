# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project( Arrangement_on_surface_2_Demo )

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL COMPONENTS Core Qt5 )

find_package( Qt5 QUIET COMPONENTS Gui Widgets)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND CGAL_Core_FOUND)
  include(${CGAL_USE_FILE})
  
  include_directories( ./ )
  # Arrangement package includes
  add_definitions(-DQT_NO_KEYWORDS)

  qt5_wrap_ui( arrangement_2_uis
    ArrangementDemoWindow.ui
    NewTabDialog.ui
    OverlayDialog.ui
    ArrangementDemoPropertiesDialog.ui
  )
  qt5_wrap_cpp( CGAL_Qt5_MOC_FILES
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
  qt5_add_resources( CGAL_Qt5_RESOURCE_FILES
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
   ${CGAL_Qt5_MOC_FILES}
   ${arrangement_2_uis}
   ${CGAL_Qt5_RESOURCE_FILES}
  )
  target_link_libraries( arrangement_2 PRIVATE
    CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core
    Qt5::Gui
  )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS arrangement_2 )

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(arrangement_2)
else()

  message(STATUS "NOTICE: This demo requires CGAL, CGAL-Core and Qt5, and will not be compiled.")

endif()

# add_subdirectory( tests )
