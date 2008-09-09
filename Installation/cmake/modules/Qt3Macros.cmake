# This file contains the following macros:
#  macro QT3_ADD_RESOURCE(outfiles inputfile ... )
#  macro QT3_AUTOMOC(inputfile ... )
#  macro QT3_GENERATE_MOC(inputfile outputfile )
#
# Adapted to Qt3 and CGAL from FindQt4.cmake (included in CMake 2.4)

INCLUDE(AddFileDependencies)


# get include dirs
MACRO (QT3_GET_MOC_INC_DIRS _moc_INC_DIRS)
    SET(${_moc_INC_DIRS})
    GET_DIRECTORY_PROPERTY(_inc_DIRS INCLUDE_DIRECTORIES)

    FOREACH(_current ${_inc_DIRS})
        SET(${_moc_INC_DIRS} "${${_moc_INC_DIRS}} -I ${_current}")
    ENDFOREACH(_current ${_inc_DIRS})
ENDMACRO(QT3_GET_MOC_INC_DIRS)


# add rule to generate ${outfile} .moc file from ${infile} (.cpp or .h)
MACRO (QT3_GENERATE_MOC infile outfile)
#     QT3_GET_MOC_INC_DIRS(moc_includes)

    GET_FILENAME_COMPONENT(infile ${infile} ABSOLUTE)
    GET_FILENAME_COMPONENT(outfile ${outfile} ABSOLUTE)

    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
        COMMAND ${QT3_MOC_EXECUTABLE}
#         ARGS ${moc_includes} -o ${outfile} ${infile}
        ARGS -o ${outfile} ${infile}
        DEPENDS ${infile})

#     ADD_FILE_DEPENDENCIES(${infile} ${outfile})
ENDMACRO (QT3_GENERATE_MOC)


# # QT3_WRAP_CPP(outfiles inputfile ... )
# MACRO (QT3_WRAP_CPP outfiles )
#     # get include dirs
#     QT3_GET_MOC_INC_DIRS(moc_includes)
#
#     FOREACH (it ${ARGN})
#         GET_FILENAME_COMPONENT(it ${it} ABSOLUTE)
#         GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
#
#         SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/moc_${outfile}.cxx)
#         ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
#             COMMAND ${QT3_MOC_EXECUTABLE}
#             ARGS ${moc_includes} -o ${outfile} ${it}
#             DEPENDS ${it})
#         SET(${outfiles} ${${outfiles}} ${outfile})
#     ENDFOREACH(it)
# ENDMACRO (QT3_WRAP_CPP)


# # QT3_WRAP_UI(outfiles inputfile ... )
# MACRO (QT3_WRAP_UI outfiles )
#
#     FOREACH (it ${ARGN})
#     GET_FILENAME_COMPONENT(outfile ${it} NAME_WE)
#     GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
#     SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/ui_${outfile}.h)
#     ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
#         COMMAND ${QT3_UIC_EXECUTABLE}
#         ARGS -o ${outfile} ${infile}
#         MAIN_DEPENDENCY ${infile})
#     SET(${outfiles} ${${outfiles}} ${outfile})
#     ENDFOREACH (it)
#
# ENDMACRO (QT3_WRAP_UI)


# QT3_ADD_RESOURCE(outfiles inputfile ... )
MACRO (QT3_ADD_RESOURCES outfiles )

    FOREACH (it ${ARGN})
    GET_FILENAME_COMPONENT(outfilename ${it} NAME_WE)
    GET_FILENAME_COMPONENT(infile ${it} ABSOLUTE)
    GET_FILENAME_COMPONENT(rc_path ${infile} PATH)
    SET(outfile ${CMAKE_CURRENT_BINARY_DIR}/qrc_${outfilename}.cxx)
    #  parse file for dependencies
    FILE(READ "${infile}" _RC_FILE_CONTENTS)
    STRING(REGEX MATCHALL "<file>[^<]*" _RC_FILES "${_RC_FILE_CONTENTS}")
    SET(_RC_DEPENDS)
    FOREACH(_RC_FILE ${_RC_FILES})
        STRING(REGEX REPLACE "^<file>" "" _RC_FILE "${_RC_FILE}")
        SET(_RC_DEPENDS ${_RC_DEPENDS} "${rc_path}/${_RC_FILE}")
    ENDFOREACH(_RC_FILE)
    ADD_CUSTOM_COMMAND(OUTPUT ${outfile}
        COMMAND ${QT3_RCC_EXECUTABLE}
        ARGS -name ${outfilename} -o ${outfile} ${infile}
        MAIN_DEPENDENCY ${infile}
        DEPENDS ${_RC_DEPENDS})
    SET(${outfiles} ${${outfiles}} ${outfile})
    ENDFOREACH (it)

ENDMACRO (QT3_ADD_RESOURCES)


# QT3_AUTOMOC(file_cpp_1 ... file_cpp_N)
#    Call this if you want to have automatic moc file handling.
#    This means if you include "foo.moc" in the source file foo.cpp
#    a moc file for the header foo.h will be created automatically.
#    if foo.h doesn't exit, the moc is created from foo.cpp
#    You can set the property SKIP_AUTOMOC using SET_SOURCE_FILES_PROPERTIES()
#    to exclude some files in the list from being processed.
MACRO(QT3_AUTOMOC)
    GET_DIRECTORY_PROPERTY(_inc_DIRS INCLUDE_DIRECTORIES)

    # For each parameter _current_FILE
    FOREACH (_current_FILE ${ARGN})
        # Get _current_FILE's full path
        GET_FILENAME_COMPONENT(_current_abs_FILE ${_current_FILE} ABSOLUTE)
        GET_FILENAME_COMPONENT(_current_abs_PATH ${_current_abs_FILE} PATH)

        # if "SKIP_AUTOMOC" is set to true, we will not handle this file here.
        GET_SOURCE_FILE_PROPERTY(_skip ${_current_abs_FILE} SKIP_AUTOMOC)
        IF ( NOT _skip AND EXISTS ${_current_abs_FILE} )
            # Read file
            FILE(READ ${_current_abs_FILE} _contents)

            STRING(REGEX MATCHALL "#include +[^ ]+\\.moc[\">]" _match "${_contents}")
            IF(_match)
                # For each "#include *.moc"
                FOREACH (_current_MOC_INCLUDE ${_match})
                    # Get name of .moc to create (full path)
                    STRING(REGEX MATCH "[^ <\"]+\\.moc" _current_MOC "${_current_MOC_INCLUDE}")
#                     SET(_moc ${_current_abs_PATH}/${_current_MOC})
                   SET(_moc ${CMAKE_CURRENT_BINARY_DIR}/${_current_MOC})

                    # Find .moc's source header (full path). The result is cached
                    # as "${_basename}_h" variable (advanced).
                    # TODO: search among headers included by ${_current_FILE}.
                    # TODO: write a macro find_file_no_cache() which does the same as
                    # find_file() without caching the result (disturbing for user).
                    GET_FILENAME_COMPONENT(_basename ${_current_MOC} NAME_WE)
#                    SET(_header ${_abs_PATH}/${_basename}.h)
                    set(_header "${_basename}_h")
                    find_file(${_header}
                              NAMES ${_basename}.h
                              PATHS ${_current_abs_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${_inc_DIRS}
                              NO_DEFAULT_PATH)
                    set ( ${_header} ${${_header}} CACHE INTERNAL "hide this" FORCE )

                    if (NOT ${_header})
                      set( moc_source "${_current_abs_FILE}" )
                    else(NOT ${_header})
                      set( moc_source ${${_header}} )
                    endif(NOT ${_header})
                    
                        # Add make rule to create .moc
#                         MESSAGE(STATUS "QT3_AUTOMOC: add rule ${_moc} <- ${moc_source}") # debug
                     include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})
                     ADD_CUSTOM_COMMAND(OUTPUT ${_moc}
                                        COMMAND ${QT3_MOC_EXECUTABLE}
#                                          ARGS ${_moc_INCS} ${_header} -o ${_moc}
                                        ${moc_source} -o ${_moc}
                                        DEPENDS ${moc_source}
                                        )
                     ADD_FILE_DEPENDENCIES(${_current_abs_FILE} ${_moc})

                ENDFOREACH (_current_MOC_INCLUDE)
            ENDIF(_match)
        ENDIF ( NOT _skip AND EXISTS ${_current_abs_FILE} )
    ENDFOREACH (_current_FILE)
ENDMACRO(QT3_AUTOMOC)

