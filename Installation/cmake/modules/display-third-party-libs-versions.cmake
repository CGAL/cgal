set(LIBRARIES_TO_CHECK
  Boost
  Ceres
  Eigen3
  GLPK
  GMP
  ITK
  ITT
  LASLIB
  libpointmatcher
  METIS
  MPFI
  MPFR
  OpenCV
  OpenGR
  OpenMesh
  OSQP
  Qt6
  SCIP
  SuiteSparse
  TBB
  Threads
  VTK
  ZLIB
)

include(${CMAKE_CURRENT_LIST_DIR}/CGAL_TweakFindBoost.cmake)

function(get_library_version header_path major_macro minor_macro patchlevel_macro version_var)
  if(EXISTS ${header_path})
    file(READ ${header_path} HEADER_CONTENT)
    if("${HEADER_CONTENT}" MATCHES "#define[ \t]+${major_macro}[ \t]+([0-9]+)")
      set(VERSION_MAJOR ${CMAKE_MATCH_1})
    endif()
    if("${HEADER_CONTENT}" MATCHES "#define[ \t]+${minor_macro}[ \t]+([0-9]+)")
      set(VERSION_MINOR ${CMAKE_MATCH_1})
    endif()
    if("${HEADER_CONTENT}" MATCHES "#define[ \t]+${patchlevel_macro}[ \t]+([0-9]+)")
      set(VERSION_PATCHLEVEL ${CMAKE_MATCH_1})
    else()
      set(VERSION_PATCHLEVEL "")
    endif()
    if(VERSION_PATCHLEVEL)
      set(${version_var} "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCHLEVEL}" PARENT_SCOPE)
    elseif(VERSION_MAJOR GREATER_EQUAL 0 OR VERSION_MINOR GREATER_EQUAL 0)
      set(${version_var} "${VERSION_MAJOR}.${VERSION_MINOR}" PARENT_SCOPE)
    endif()
  endif()
endfunction()

function(check_library cgal_3rdparty_lib)
  set(QT_NO_CREATE_VERSIONLESS_TARGETS ON)
  set(CMAKE_FIND_PACKAGE_QUIET TRUE)
  string(TOUPPER ${cgal_3rdparty_lib} cgal_3rdparty_lib_upper)
  find_package(${cgal_3rdparty_lib} QUIET)
  set(CMAKE_FIND_PACKAGE_QUIET FALSE)
  if(${cgal_3rdparty_lib}_FOUND)
    set(version_var "")
    if(DEFINED ${cgal_3rdparty_lib}_VERSION)
      set(version_var ${${cgal_3rdparty_lib}_VERSION})
    elseif(DEFINED ${cgal_3rdparty_lib}_VERSION_STRING)
      set(version_var ${${cgal_3rdparty_lib}_VERSION_STRING})
    elseif(DEFINED ${cgal_3rdparty_lib_upper}_VERSION)
      set(version_var ${${cgal_3rdparty_lib_upper}_VERSION})
    elseif(DEFINED ${cgal_3rdparty_lib_upper}_VERSION_STRING)
      set(version_var ${${cgal_3rdparty_lib_upper}_VERSION_STRING})
    elseif(${cgal_3rdparty_lib} STREQUAL "GMP")
      set(version_var "")
      get_library_version("${GMP_INCLUDE_DIR}/gmp.h" "__GNU_MP_VERSION" "__GNU_MP_VERSION_MINOR" "__GNU_MP_VERSION_PATCHLEVEL" version_var)
      if(NOT version_var)
        file(READ "${GMP_INCLUDE_DIR}/gmp.h" GMP_HEADER_CONTENT)
        string(REGEX MATCHALL "#include[ \t]+\"([^\"]+)\"" INCLUDED_HEADERS "${GMP_HEADER_CONTENT}")
        foreach(INCLUDED_HEADER ${INCLUDED_HEADERS})
          string(REGEX REPLACE "#include[ \t]+\"([^\"]+)\"" "\\1" GMP_ARCH_HEADER "${INCLUDED_HEADER}")
          set(GMP_ARCH_HEADER_PATH "${GMP_INCLUDE_DIR}/${GMP_ARCH_HEADER}")
          if(EXISTS ${GMP_ARCH_HEADER_PATH})
            get_library_version("${GMP_ARCH_HEADER_PATH}" "__GNU_MP_VERSION" "__GNU_MP_VERSION_MINOR" "__GNU_MP_VERSION_PATCHLEVEL" version_var)
            if(version_var)
              break()
            endif()
          endif()
        endforeach()
      endif()
    elseif(${cgal_3rdparty_lib} STREQUAL "GLPK")
      get_library_version("${GLPK_INCLUDE_DIR}/glpk.h" "GLP_MAJOR_VERSION" "GLP_MINOR_VERSION" "" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "MPFR")
      get_library_version("${MPFR_INCLUDE_DIR}/mpfr.h" "MPFR_VERSION_MAJOR" "MPFR_VERSION_MINOR" "MPFR_VERSION_PATCHLEVEL" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "METIS")
      get_library_version("${METIS_INCLUDE_DIR}/metis.h" "METIS_VER_MAJOR" "METIS_VER_MINOR" "METIS_VER_SUBMINOR" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "SuiteSparse")
      get_library_version("${SuiteSparse_Config_INCLUDE_DIR}/SuiteSparse_config.h" "SUITESPARSE_MAIN_VERSION" "SUITESPARSE_SUB_VERSION" "SUITESPARSE_SUBSUB_VERSION" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "LASLIB")
      get_library_version("${LASLIB_INCLUDE_DIR}/lasdefinitions.hpp" "LAS_TOOLS_VERSION" "" "" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "ITT")
      get_library_version("${ITT_INCLUDE_DIR}/ittnotify.h" "ITT_MAJOR" "ITT_MINOR" "" version_var)
    elseif(${cgal_3rdparty_lib} STREQUAL "OpenMesh")
      if (TARGET OpenMeshCore)
        get_target_property(OpenMesh_INCLUDE_DIRS OpenMeshCore INTERFACE_INCLUDE_DIRECTORIES)
        set(CONFIG_FILE_PATH "${OpenMesh_INCLUDE_DIRS}/OpenMesh/Core/System/config.h")
        if(EXISTS ${CONFIG_FILE_PATH})
          file(READ ${CONFIG_FILE_PATH} CONFIG_CONTENT)
          if("${CONFIG_CONTENT}" MATCHES "#define[ \t]+OM_VERSION[ \t]+(0x[0-9A-F]+)")
            set(VERSION_HEX ${CMAKE_MATCH_1})
            math(EXPR VERSION_MAJOR "(${VERSION_HEX} & 0xF0000) >> 16")
            math(EXPR VERSION_MINOR "(${VERSION_HEX} & 0x0FF00) >> 8")
            math(EXPR VERSION_PATCH "(${VERSION_HEX} & 0x000FF)")
            set(version_var "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
          endif()
        endif()
      endif()
    endif()
    if(version_var)
      message(STATUS "Third-party library ${cgal_3rdparty_lib} ${version_var}")
    else()
      message(STATUS "Third-party library ${cgal_3rdparty_lib} found")
    endif()
  else()
    message(STATUS "Third-party library ${cgal_3rdparty_lib} not found")
  endif()
endfunction()

foreach(cgal_3rdparty_lib IN LISTS LIBRARIES_TO_CHECK)
  check_library(${cgal_3rdparty_lib})
endforeach()
