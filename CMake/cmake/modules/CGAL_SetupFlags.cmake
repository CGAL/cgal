set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CGAL_CXX_FLAGS}")  
message( STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}" )

if ( CGAL_RELEASE )
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CGAL_CXX_FLAGS_RELEASE}")  
  message( STATUS "Release compiler flags: ${CMAKE_CXX_FLAGS_RELEASE}" )
endif()

if ( CGAL_DEBUG )
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CGAL_CXX_FLAGS_DEBUG}") 
  message( STATUS "Debug compiler flags: ${CMAKE_CXX_FLAGS_DEBUG}" )
endif()

if ( BUILD_SHARED_LIBS )

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${CGAL_SHARED_LINKER_FLAGS}")
  message( STATUS "Shared linker flags: ${CMAKE_SHARED_LINKER_FLAGS}" )
  
  if ( CGAL_RELEASE )
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} ${CGAL_SHARED_LINKER_FLAGS_RELEASE}")
    message( STATUS "Release shared linker flags: ${CMAKE_SHARED_LINKER_FLAGS_RELEASE}" )
  endif()
  
  if ( CGAL_DEBUG )
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${CGAL_SHARED_LINKER_FLAGS_DEBUG}")
    message( STATUS "Debug shared linker flags: ${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
  endif()
  
else()

  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CGAL_MODULE_LINKER_FLAGS}")
  message( STATUS "Module linker flags: ${CMAKE_MODULE_LINKER_FLAGS}" )
  
  if ( CGAL_RELEASE )
    set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS_RELEASE} ${CGAL_MODULE_LINKER_FLAGS_RELEASE}")
    message( STATUS "Release module linker flags: ${CMAKE_MODULE_LINKER_FLAGS_RELEASE}" )
  endif()
  
  if ( CGAL_DEBUG )
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "${CMAKE_MODULE_LINKER_FLAGS_DEBUG} ${CGAL_MODULE_LINKER_FLAGS_DEBUG}")
    message( STATUS "Debug module linker flags: ${CMAKE_MODULE_LINKER_FLAGS_DEBUG}" )
  endif()
  
endif()

add_subdirectory(src)

set( CGAL_LIBRARIES ${CGAL_CORE_LIBRARY} ${CGAL_LIBRARY} ${CGAL_IMAGE_IO_LIBRARY} ${CGAL_PDB_LIBRARY} ${CGAL_QT_LIBRARY} )

hide_variable(EXECUTABLE_OUTPUT_PATH)
hide_variable(LIBRARY_OUTPUT_PATH)

#--------------------------------------------------------------------------------------------------
#
#                                    -= USER SIDE SETTINGS =-
#
#--------------------------------------------------------------------------------------------------

# FindCGAL and UseCGAL are platform specific so they are generated and stored in the binary folder.
configure_file(${CGAL_SOURCE_DIR}/CGALConfig_binary.cmake.in  ${CGAL_BINARY_DIR}/CGALConfig.cmake       @ONLY IMMEDIATE)



if ( SOURCE_INSTALL )
  configure_file(${CGAL_SOURCE_DIR}/CGALConfig_install.cmake.source.in ${CGAL_BINARY_DIR}/cmake/CGALConfig.cmake @ONLY IMMEDIATE)
else()
  configure_file(${CGAL_SOURCE_DIR}/CGALConfig_install.cmake.fhs.in    ${CGAL_BINARY_DIR}/cmake/CGALConfig.cmake @ONLY IMMEDIATE)
endif()


#--------------------------------------------------------------------------------------------------
#
#                                    -= Installation Commands =-
#
#--------------------------------------------------------------------------------------------------

if ( WITH_INSTALL )

  # WARNING: Use only relative paths; full paths break CPack!
  # DESTINATION option is mandatory; skipping it breaks CPack!

  install(FILES CHANGES "INSTALL" INSTALL.MacOSX INSTALL.win32.txt LICENSE LICENSE.FREE_USE LICENSE.LGPL LICENSE.QPL README "VERSION"  
          DESTINATION ${CGAL_INSTALL_DOC_DIR}
         )


  install(DIRECTORY include/CGAL                      DESTINATION ${CGAL_INSTALL_INC_DIR}    )
  install(DIRECTORY "${CGAL_BINARY_DIR}/include/CGAL" DESTINATION ${CGAL_INSTALL_INC_DIR}    )

  install(DIRECTORY scripts/                          DESTINATION ${CGAL_INSTALL_BIN_DIR}    )

  install(DIRECTORY cmake/modules/                        DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )
  install(FILES cmake/modules/UseCGAL.cmake               DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )

  if ( GMP_IN_AUXILIARY )
    install(DIRECTORY auxiliary/gmp/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
    install(DIRECTORY auxiliary/gmp/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
  endif()

  if ( TAUCS_IN_AUXILIARY )
    install(DIRECTORY auxiliary/taucs/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
    install(DIRECTORY auxiliary/tacus/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
  endif()

  if ( ZLIB_IN_AUXILIARY )
    install(DIRECTORY auxiliary/zlib/include/ DESTINATION ${CGAL_INSTALL_INC_DIR} )
    install(DIRECTORY auxiliary/zlib/lib/     DESTINATION ${CGAL_INSTALL_LIB_DIR} )
  endif()


  if ( SOURCE_INSTALL )
    install(FILES ${CGAL_BINARY_DIR}/cmake/CGALConfig.cmake DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
    install(FILES     CMakeLists.txt                        DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
    install(DIRECTORY auxiliary                             DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
    install(DIRECTORY src                                   DESTINATION ${CGAL_INSTALL_ROOT_DIR} )
  else()
    install(FILES ${CGAL_BINARY_DIR}/cmake/CGALConfig.cmake DESTINATION ${CGAL_INSTALL_CMAKE_DIR} )
  endif()
  
endif()


#--------------------------------------------------------------------------------------------------
#
#                                    -= APPLICATIONS =-
#
#--------------------------------------------------------------------------------------------------


# This allows programs to locate CGALConfig.cmake
set(CGAL_DIR ${CGAL_BINARY_DIR} )

# This script does not configure demos and examples by default, but its does if requested.
optional_add_subdirectory( demo     ON EXCLUDE_FROM_ALL )
optional_add_subdirectory( examples ON EXCLUDE_FROM_ALL )


#--------------------------------------------------------------------------------------------------
#
#                                            -= CPack =-
#
#--------------------------------------------------------------------------------------------------

option( WITH_CPACK "Create package generation rules")
if( WITH_CPACK AND EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake" )

    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CGAL - Computational Geometry Algorithms Library")
    set(CPACK_PACKAGE_VENDOR "CGAL Open Source Project")
    set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    set(CPACK_PACKAGE_VERSION_MAJOR "${CGAL_MAJOR_VERSION}")
    set(CPACK_PACKAGE_VERSION_MINOR "${CGAL_MINOR_VERSION}")
    set(CPACK_PACKAGE_VERSION_PATCH "${CGAL_BUILD_VERSION}")
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "CGAL ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "CGAL-${CGAL_VERSION}")
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
    
    if(NOT DEFINED CPACK_SYSTEM_NAME)
      set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
    endif()
    
    if(${CPACK_SYSTEM_NAME} MATCHES Windows)
      if(CMAKE_CL_64)
        set(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
      else()
        set(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
      endif()
    endif()
    
    if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
      set(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
    endif()
    
    set(CPACK_PACKAGE_EXECUTABLES "CGAL" "CGAL")

    if(WIN32 AND NOT UNIX)
        set(CPACK_GENERATOR "NSIS")
        # There is a bug in NSI that does not handle full unix paths properly. Make
        # sure there is at least one set of four (4) backlasshes.
        #set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}\\\\cgal_install.gif")
        #set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\CGAL.exe")
        set(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} Computational Geometry Algorithms Library")
        set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.cgal.org")
        set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.cgal.com")
        set(CPACK_NSIS_CONTACT "info@cgal.org")
        set(CPACK_NSIS_MODifY_PATH ON)
    else()
        set(CPACK_STRIP_FILES "")
        set(CPACK_SOURCE_STRIP_FILES "")
    endif()

    INCLUDE(CPack)
    
endif()

