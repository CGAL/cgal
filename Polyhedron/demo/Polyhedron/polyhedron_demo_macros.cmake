include(AddFileDependencies)

  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      if(NOT ${option} STREQUAL "NO_MOC")
        set(other_sources ${ARGN})
        set(option "")
      else()
        set(other_sources ${ARGN_TAIL})
      endif()
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    if("${option}" STREQUAL "NO_MOC") 
      set(option "")
      set(moc_file_name "")
    else()
      set(moc_file_name ${plugin_implementation_base_name}.moc )
      qt4_generate_moc( ${plugin_implementation_base_name}.cpp "${CMAKE_CURRENT_BINARY_DIR}/${moc_file_name}" )
      add_file_dependencies( ${moc_file_name} "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )
    endif()

    add_library(${plugin_name} MODULE ${option} ${moc_file_name} ${plugin_implementation_base_name}.cpp ${other_sources})
    set_property(TARGET ${plugin_name}
      PROPERTY LIBRARY_OUTPUT_DIRECTORY
      "${CGAL_POLYHEDRON_DEMO_PLUGINS_DIR}")

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with the demo_framework
    if(TARGET demo_framework)
      target_link_libraries( ${plugin_name} demo_framework)
    else()
      target_link_libraries( ${plugin_name} Polyhedron_demo_framework)
    endif()
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
    add_dependencies( ${plugin_name} Polyhedron_3 )
  endmacro(polyhedron_demo_plugin)
