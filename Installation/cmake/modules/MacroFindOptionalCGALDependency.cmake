include(MacroOptionalFindPackage)

macro (find_optional_cgal_dependency _name )
  macro_optional_find_package(${_name})
  if(WITH_${_name} AND ${_name}_FOUND )
    set(CGAL_3RD_PARTY_INCLUDE_DIRS    ${CGAL_3RD_PARTY_INCLUDE_DIRS}    ${${_name}_INCLUDE_DIR}   CACHE INTERNAL "" FORCE )
    set(CGAL_3RD_PARTY_LIBRARIES_DIRS  ${CGAL_3RD_PARTY_LIBRARIES_DIRS}  ${${_name}_LIBRARIES_DIR} CACHE INTERNAL "" FORCE )
    set(CGAL_3RD_PARTY_LIBRARIES       ${CGAL_3RD_PARTY_LIBRARIES}       ${${_name}_LIBRARIES}     CACHE INTERNAL "" FORCE )
    set(CGAL_USE_${_name} 1)
  endif()
endmacro()

