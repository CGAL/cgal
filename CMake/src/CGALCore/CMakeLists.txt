# src/Core/CMakeLists.txt

# CORE needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)
    # Make sure the compiler can find GMP include files
    include_directories (${GMP_INCLUDE_DIR})

    # Make sure the compiler can find CORE include files
    include_directories (BEFORE ${CGAL_SOURCE_DIR}/include)

    # Folder where the compiler will put the CORE library
    set(LIBRARY_OUTPUT_PATH ${CGAL_BINARY_DIR}/lib)

    # Get list of source files in this folder
    AUX_SOURCE_DIRECTORY(. SOURCES)

    # Create core++ library
    ADD_LIBRARY(core++ ${SOURCES})

    # Link it with GMP
    target_link_libraries(core++ ${GMP_LIBRARIES})

    # Install library
    # WARNING: Use only relative paths; full paths break CPack!
    #          DESTINATION option is mandatory; skipping it breaks CPack!
    INSTALL(TARGETS core++
	    RUNTIME DESTINATION ${CGAL_BIN_INSTALL_DIR}
	    LIBRARY DESTINATION ${CGAL_LIB_INSTALL_DIR}
	    ARCHIVE DESTINATION ${CGAL_LIB_INSTALL_DIR})
endif(GMP_FOUND)
