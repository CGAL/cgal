# src/CGAL/CMakeLists.txt
#
# Build CGAL library.
#
# Make sure the compiler can find include files
include_directories (${CGAL_INCLUDE_DIRS})

# Add compiler flags CGAL needs to work on this platform.
ADD_DEFINITIONS(${CGAL_LIB_CXX_FLAGS})
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CGAL_LIB_EXE_LINKER_FLAGS}")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CGAL_LIB_SHARED_LINKER_FLAGS}")

# List of source files for CGAL library
set(CGAL_SOURCES
    assertions.cpp
    Bbox_2_intersections.cpp
    Bbox_3_intersections.cpp
    cgal_logo.cpp
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

# Create a library called "CGAL" from ${CGAL_SOURCES} source files.
# The extension is automatically found.
add_library (CGAL ${CGAL_SOURCES})

# Link it with third party libraries
target_link_libraries(CGAL ${CGAL_3RD_PARTY_LIBRARIES})
if (CGAL_USE_CGAL_CORE AND NOT AUTO_LINK_ENABLED )
  add_dependencies(CGAL core++)
endif()

# Apply user-defined properties to the library targets.
if(CGAL_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(CGAL PROPERTIES ${CGAL_LIBRARY_PROPERTIES})
endif()

# Install library
# WARNING: Use only relative paths; full paths break CPack!
#          DESTINATION option is mandatory; skipping it breaks CPack!
INSTALL(TARGETS CGAL
	      RUNTIME DESTINATION ${CGAL_BIN_INSTALL_DIR}
	      LIBRARY DESTINATION ${CGAL_LIB_INSTALL_DIR}
       	ARCHIVE DESTINATION ${CGAL_LIB_INSTALL_DIR}
        )

