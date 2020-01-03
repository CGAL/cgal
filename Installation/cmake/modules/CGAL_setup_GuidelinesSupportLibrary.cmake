if(NOT TARGET cxx::gsl)
  include(ExternalProject)

  if(NOT CGAL_NO_EXTERNAL_PROJECTS)
    ExternalProject_Add(gsl-download
      URL               https://github.com/microsoft/GSL/archive/v2.0.0.tar.gz
      DOWNLOAD_NAME     GSL-2.0.0.tar.gz
      URL_HASH          SHA512=7339527222c8a97a94c0bb4038b3d142045ec5d80995e628574ac96f4d9d13c41ad70fbe0d8390586dc0db8d9ea55107dbc95de80f7335eb78ef9d2e7047d726
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ""
      INSTALL_COMMAND   ""
      )
    ExternalProject_Get_Property(gsl-download source_dir)
    find_path(gsl_include_dir
      NAMES gsl/gsl
      PATHS ${source_dir}/include
      DOC "Include directory of the C++ Guideline Support Library")
  endif()
  if(NOT gsl_include_dir)
    find_path(gsl_include_dir
      NAMES gsl/gsl
      PATH_SUFFIXES guidelines-support-library
      DOC "Include directory of the C++ Guideline Support Library (version 2.0 or later)")
  endif()
  if(gsl_include_dir)
    add_library(CGAL-gsl INTERFACE)
    target_include_directories(CGAL-gsl SYSTEM INTERFACE ${gsl_include_dir})
    if(TARGET gsl-download)
      add_dependencies(CGAL-gsl gsl-download)
    endif()
    add_library(cxx::gsl ALIAS CGAL-gsl)
  endif()
endif()
