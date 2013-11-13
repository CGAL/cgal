# fill string with spaces
macro (acg_format_string str length return)
    string (LENGTH "${str}" _str_len)
    math (EXPR _add_chr "${length} - ${_str_len}")
    set (${return} "${str}")
    while (_add_chr GREATER 0)
        set (${return} "${${return}} ")
        math (EXPR _add_chr "${_add_chr} - 1")
    endwhile ()
endmacro ()

# print message with color escape sequences if CMAKE_COLOR_MAKEFILE is set
string (ASCII 27 _escape)
function (acg_color_message _str)
    if (CMAKE_COLOR_MAKEFILE AND NOT WIN32)
        message (${_str})
    else ()
        string (REGEX REPLACE "${_escape}.[0123456789;]*m" "" __str ${_str})
        message (${__str})
    endif ()
endfunction ()

# info header
function (acg_print_configure_header _id _name)
    acg_format_string ("${_name}" 40 _project)
    acg_format_string ("${${_id}_VERSION}" 40 _version)
    acg_color_message ("\n${_escape}[40;37m************************************************************${_escape}[0m")
    acg_color_message ("${_escape}[40;37m* ${_escape}[1;31mACG ${_escape}[0;40;34mBuildsystem${_escape}[0m${_escape}[40;37m                                          *${_escape}[0m")
    acg_color_message ("${_escape}[40;37m*                                                          *${_escape}[0m")
    acg_color_message ("${_escape}[40;37m* Package : ${_escape}[32m${_project} ${_escape}[37m      *${_escape}[0m")
    acg_color_message ("${_escape}[40;37m* Version : ${_escape}[32m${_version} ${_escape}[37m      *${_escape}[0m")
    acg_color_message ("${_escape}[40;37m************************************************************${_escape}[0m")
endfunction ()

# info line
function (acg_print_configure_footer)
    acg_color_message ("${_escape}[40;37m************************************************************${_escape}[0m\n")
endfunction ()
