# src/CGAL/CMakeLists.txt
#
# The following variables must be set from the parent script:
# The first two cannot be empty
#
# CGAL_BINARY_DIR
# CGAL_LIB_INSTALL_DIR
#
# CGAL_INCLUDE_DIRS
# CGAL_DEFINITIONS
# CGAL_3RD_PARTY_INCLUDE_DIRS
# CGAL_3RD_PARTY_LIBRARIES
# CGAL_3RD_PARTY_DEFINITIONS
# CGAL_LIB_SHARED_LINKER_FLAGS
# CGAL_VERSION
# CGAL_MAJOR_VERSION
# CGAL_MINOR_VERSION

   
assert(CGAL_BINARY_DIR)
assert(CGAL_LIB_INSTALL_DIR)
   
include_directories (${CGAL_INCLUDE_DIRS})
include_directories (${CGAL_3RD_PARTY_INCLUDE_DIRS})

add_definitions(${CGAL_DEFINITIONS})
add_definitions(${CGAL_3RD_PARTY_DEFINITIONS})

set(LIBRARY_OUTPUT_PATH ${CGAL_BINARY_DIR}/lib)

if (WIN32)
  set(BUILD_SHARED_LIBS OFF)
endif()

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CGAL_LIB_SHARED_LINKER_FLAGS}")

# List of source files for CGAL library
set(SOURCES
    assertions.cpp
    Bbox_2_intersections.cpp
    Bbox_3_intersections.cpp
    Color.cpp
    File_header_extended_OFF.cpp
    File_header_OFF.cpp
    File_scanner_OFF.cpp
    File_writer_inventor.cpp
    File_writer_OFF.cpp
    File_writer_VRML_2.cpp
    File_writer_wavefront.cpp
    Geomview_stream.cpp
    Interval_arithmetic.cpp
    io.cpp
    JAMA_numeric_solver.cpp
    KDS_Log.cpp
    kernel.cpp
    MP_Float.cpp
    numeric_solvers_support.cpp
    NefPolynomial.cpp
    Random.cpp
    Real_timer.cpp
    Timer.cpp
    Turkowski_numeric_solver.cpp
)

add_library (CGAL ${SOURCES})

target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})

if (CGAL_USE_CGAL_CORE)
  add_dependencies(CGAL CGAL_CORE)
endif()

set_target_properties(CGAL PROPERTIES VERSION "${CGAL_VERSION}" SOVERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")

if ( AUTO_LINK_ENABLED )

  set(CGAL_RELEASE_MANGLED_NAME        "cgal${TOOLSET}-mt" ) 
  set(CGAL_DEBUG_MANGLED_NAME          "cgal${TOOLSET}-mt-gd" ) 
  set(CGAL_MINSIZEREL_MANGLED_NAME     "cgal${TOOLSET}-mt-o" ) 
  set(CGAL_RELWITHDEBINFO_MANGLED_NAME "cgal${TOOLSET}-mt-g" ) 
  
  set_target_properties( CGAL PROPERTIES DEBUG_OUTPUT_NAME          "${CGAL_DEBUG_MANGLED_NAME}" 
                                         RELEASE_OUTPUT_NAME        "${CGAL_RELEASE_MANGLED_NAME}" 
                                         MINSIZEREL_OUTPUT_NAME     "${CGAL_MINSIZEREL_MANGLED_NAME}" 
                                         RELWITHDEBINFO_OUTPUT_NAME "${CGAL_RELWITHDEBINFO_MANGLED_NAME}" 
                       )
  
  add_custom_command(TARGET CGAL POST_BUILD COMMAND move /Y $(TargetDir)$(TargetName)$(TargetExt) ../../lib)
  add_custom_command(TARGET CGAL POST_BUILD COMMAND if exist "$(TargetDir)$(TargetName).pdb." move /Y $(TargetDir)$(TargetName).pdb ../../lib)
  
endif()

install(TARGETS CGAL DESTINATION ${CGAL_LIB_INSTALL_DIR} )

