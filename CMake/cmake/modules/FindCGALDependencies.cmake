include(MacroFindOptionalCGALDependency)

find_package(Boost REQUIRED)
if(Boost_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS}   ${Boost_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${Boost_LIBRARY_DIR})
    set(CGAL_USE_BOOST 1)
endif()

find_package(BOOST_THREAD REQUIRED)
if(BOOST_THREAD_FOUND)
    if ( NOT AUTO_LINK_ENABLED )
        set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${BOOST_THREAD_LIBRARIES})
    endif( NOT AUTO_LINK_ENABLED )
    set(CGAL_USE_BOOST_THREAD 1)
endif()

find_package(OpenGL)
if(OPENGL_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})
    set(CGAL_3RD_PARTY_LIBRARIES    ${CGAL_3RD_PARTY_LIBRARIES}    ${OPENGL_LIBRARIES})
    set(CGAL_USE_OPENGL 1)
endif()

#
# find_optional_cgal_dependency(ABC) uses the option WITH_ABC to select or skip the dependency.
# If found:
#   ABC_FOUND is set to TRUE
#   CGAL_USE_ABC is set to 1
#   CGAL_3RD_PARTY_INCLUDE_DIRS    is added with ABC_INCLUDE_DIR
#   CGAL_3RD_PARTY_LIBRARIES_DIRS  is added with ABC_LIBRARIES_DIR
#   CGAL_3RD_PARTY_LIBRARIES       is added with ABC_LIBRARIES

find_optional_cgal_dependency(GMP)
find_optional_cgal_dependency(MPFR)
find_optional_cgal_dependency(ZLIB)
find_optional_cgal_dependency(TAUCS)

if ( NOT WIN32 )
  find_optional_cgal_dependency(GMPXX)
endif()


macro_optional_find_package(CGAL_CORE)
if(WITH_CGAL_CORE AND CGAL_CORE_FOUND )
  set(CGAL_USE_CGAL_CORE 1)
endif()

find_package(Qt3)

get_dependency_version(BOOST Boost )
get_dependency_version(OPENGL OpenGL)
get_dependency_version(GMP)
get_dependency_version(MPFR)
get_dependency_version(ZLIB)
get_dependency_version(TAUCS)
get_dependency_version(QT)

