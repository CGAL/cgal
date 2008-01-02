include(MacroFindOptionalCGALDependency)

find_package(Boost2 REQUIRED)
if(Boost_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS}   ${Boost_INCLUDE_DIRS})
    set(CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${Boost_LIBRARY_DIRS})
    set(CGAL_USE_BOOST 1)
endif()

find_package(BOOST_THREAD REQUIRED)
if(BOOST_THREAD_FOUND)
    set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS}   ${Boost_INCLUDE_DIRS})
    set(CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${Boost_LIBRARY_DIRS})
    set(CGAL_3RD_PARTY_LIBRARIES      ${CGAL_3RD_PARTY_LIBRARIES}      ${BOOST_THREAD_LIBRARIES})
    set(CGAL_USE_BOOST_THREAD 1)
message( STATUS "BOOST_THREAD_LIBRARIES=${BOOST_THREAD_LIBRARIES}" )    
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
find_optional_cgal_dependency(GMPXX)
find_optional_cgal_dependency(CGAL_CORE)

if ( NOT CGAL_CORE_FOUND )
  find_optional_cgal_dependency(CORE)
endif()

# TODO: Write FindLIDIA.cmake
# TODO: Write FindLEDA.cmake
# TODO: Write FindLEDAWIN.cmake

