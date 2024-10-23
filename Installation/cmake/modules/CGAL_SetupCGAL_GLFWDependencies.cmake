if(CGAL_SetupCGAL_GLFWDependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_GLFWDependencies_included TRUE)

include(CMakeDependentOption)

if (UNIX AND NOT APPLE)
  if(DEFINED ENV{XDG_SESSION_TYPE})
    set(XDG_SESSION_TYPE "$ENV{XDG_SESSION_TYPE}")
  else()
    set(XDG_SESSION_TYPE "")
  endif()
  if (XDG_SESSION_TYPE STREQUAL "wayland")
    set(GLFW_USE_WAYLAND ON)
  else()
    set(GLFW_USE_WAYLAND OFF)
  endif()
endif()

if (GLFW_USE_WAYLAND)
  add_definitions(-DGLFW_USE_WAYLAND)
endif()

# --------------------------------------------------
# From glfw/CMakeLists.txt and glfw/src/CMakeLists
# --------------------------------------------------

cmake_dependent_option(GLFW_BUILD_WIN32 "Build support for Win32" ON "WIN32" OFF)
cmake_dependent_option(GLFW_BUILD_COCOA "Build support for Cocoa" ON "APPLE" OFF)
cmake_dependent_option(GLFW_BUILD_X11 "Build support for X11" ON "UNIX; NOT APPLE; NOT GLFW_USE_WAYLAND" OFF)
cmake_dependent_option(GLFW_BUILD_WAYLAND "Build support for Wayland" ON "UNIX; NOT APPLE; GLFW_USE_WAYLAND" OFF)

cmake_dependent_option(GLFW_USE_HYBRID_HPG "Force use of high-performance GPU on hybrid systems" OFF
                       "WIN32" OFF)
cmake_dependent_option(USE_MSVC_RUNTIME_LIBRARY_DLL "Use MSVC runtime library DLL" ON
                       "MSVC" OFF)

list(APPEND CMAKE_MODULE_PATH .)

find_package(Threads REQUIRED)

#--------------------------------------------------------------------
# Report backend selection
#--------------------------------------------------------------------
if (GLFW_BUILD_WIN32)
    message(STATUS "Including Win32 support")
endif()
if (GLFW_BUILD_COCOA)
    message(STATUS "Including Cocoa support")
endif()
if (GLFW_BUILD_WAYLAND)
    message(STATUS "Including Wayland support")
endif()
if (GLFW_BUILD_X11)
    message(STATUS "Including X11 support")
endif()

#--------------------------------------------------------------------
# Apply Microsoft C runtime library option
#--------------------------------------------------------------------
if (MSVC AND NOT USE_MSVC_RUNTIME_LIBRARY_DLL)
  if (CMAKE_VERSION VERSION_LESS 3.15)
    foreach (flag CMAKE_C_FLAGS
                  CMAKE_C_FLAGS_DEBUG
                  CMAKE_C_FLAGS_RELEASE
                  CMAKE_C_FLAGS_MINSIZEREL
                  CMAKE_C_FLAGS_RELWITHDEBINFO)

      if (flag MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag} "${${flag}}")
      endif()
      if (flag MATCHES "/MDd")
        string(REGEX REPLACE "/MDd" "/MTd" ${flag} "${${flag}}")
      endif()

    endforeach()
  else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
endif()


# Define directories for each dependency
set(GLFW_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glfw")
set(GLAD_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glad")
set(STB_IMAGE_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glfw/deps")

# Include directories for compilation
include_directories(
  ${GLFW_SOURCE_DIR}/include # for glfw3.h
  ${GLAD_SOURCE_DIR}/include # for glad.h
  ${STB_IMAGE_SOURCE_DIR}    # for stb_image_write.h 
)

file(GLOB GLFW_SOURCES ${GLFW_SOURCE_DIR}/src/internal.h
                       ${GLFW_SOURCE_DIR}/src/platform.h
                       ${GLFW_SOURCE_DIR}/src/mappings.h
                       ${GLFW_SOURCE_DIR}/src/context.c
                       ${GLFW_SOURCE_DIR}/src/init.c
                       ${GLFW_SOURCE_DIR}/src/input.c
                       ${GLFW_SOURCE_DIR}/src/monitor.c
                       ${GLFW_SOURCE_DIR}/src/platform.c
                       ${GLFW_SOURCE_DIR}/src/vulkan.c
                       ${GLFW_SOURCE_DIR}/src/window.c
                       ${GLFW_SOURCE_DIR}/src/egl_context.c
                       ${GLFW_SOURCE_DIR}/src/osmesa_context.c
                       ${GLFW_SOURCE_DIR}/src/null_platform.h
                       ${GLFW_SOURCE_DIR}/src/null_joystick.h
                       ${GLFW_SOURCE_DIR}/src/null_init.c
                       ${GLFW_SOURCE_DIR}/src/null_monitor.c
                       ${GLFW_SOURCE_DIR}/src/null_window.c
                       ${GLFW_SOURCE_DIR}/src/null_joystick.c
)
add_library(glfw ${GLFW_SOURCES})

if (APPLE)
  file(GLOB GLFW_APPLE_SOURCES ${GLFW_SOURCE_DIR}/src/cocoa_time.h
                               ${GLFW_SOURCE_DIR}/src/cocoa_time.c
                               ${GLFW_SOURCE_DIR}/src/posix_thread.h
                               ${GLFW_SOURCE_DIR}/src/posix_thread.c
                               ${GLFW_SOURCE_DIR}/src/posix_module.c
  )
  target_sources(glfw PRIVATE ${GLFW_APPLE_SOURCES})
elseif (WIN32)
  file(GLOB GLFW_WIN32_SOURCES ${GLFW_SOURCE_DIR}/src/win32_time.h
                               ${GLFW_SOURCE_DIR}/src/win32_time.c
                               ${GLFW_SOURCE_DIR}/src/win32_thread.h
                               ${GLFW_SOURCE_DIR}/src/win32_thread.c
                               ${GLFW_SOURCE_DIR}/src/win32_module.c
  )
  target_sources(glfw PRIVATE ${GLFW_WIN32_SOURCES})
else()
  file(GLOB GLFW_UNIX_SOURCES ${GLFW_SOURCE_DIR}/src/posix_time.h
                               ${GLFW_SOURCE_DIR}/src/posix_time.c
                               ${GLFW_SOURCE_DIR}/src/posix_thread.h
                               ${GLFW_SOURCE_DIR}/src/posix_thread.c
                               ${GLFW_SOURCE_DIR}/src/posix_module.c
  )
  target_sources(glfw PRIVATE ${GLFW_UNIX_SOURCES})
endif()

if (GLFW_BUILD_COCOA)
  file(GLOB GLFW_COCOA_SOURCES ${GLFW_SOURCE_DIR}/src/cocoa_platform.h
                               ${GLFW_SOURCE_DIR}/src/cocoa_joystick.h
                               ${GLFW_SOURCE_DIR}/src/cocoa_joystick.m
                               ${GLFW_SOURCE_DIR}/src/cocoa_init.m
                               ${GLFW_SOURCE_DIR}/src/cocoa_monitor.m
                               ${GLFW_SOURCE_DIR}/src/cocoa_window.m
                               ${GLFW_SOURCE_DIR}/src/nsgl_context.m
  )
  target_compile_definitions(glfw PRIVATE _GLFW_COCOA)
  target_sources(glfw PRIVATE ${GLFW_COCOA_SOURCES})
endif()

if (GLFW_BUILD_WIN32)
  file(GLOB GLFW_WIN32_SOURCES ${GLFW_SOURCE_DIR}/src/win32_platform.h
                               ${GLFW_SOURCE_DIR}/src/win32_joystick.h
                               ${GLFW_SOURCE_DIR}/src/win32_joystick.c
                               ${GLFW_SOURCE_DIR}/src/win32_init.c
                               ${GLFW_SOURCE_DIR}/src/win32_monitor.c
                               ${GLFW_SOURCE_DIR}/src/win32_window.c
                               ${GLFW_SOURCE_DIR}/src/wgl_context.c
  )
  target_compile_definitions(glfw PRIVATE _GLFW_WIN32)
  target_sources(glfw PRIVATE ${GLFW_WIN32_SOURCES})
endif()

if (GLFW_BUILD_X11)
  file(GLOB GLFW_X11_SOURCES ${GLFW_SOURCE_DIR}/src/x11_platform.h
                             ${GLFW_SOURCE_DIR}/src/xkb_unicode.h
                             ${GLFW_SOURCE_DIR}/src/xkb_unicode.c
                             ${GLFW_SOURCE_DIR}/src/x11_init.c
                             ${GLFW_SOURCE_DIR}/src/x11_monitor.c
                             ${GLFW_SOURCE_DIR}/src/x11_window.c
                             ${GLFW_SOURCE_DIR}/src/glx_context.c
  )
  target_compile_definitions(glfw PRIVATE _GLFW_X11)
  target_sources(glfw PRIVATE ${GLFW_X11_SOURCES})
endif()

if (GLFW_BUILD_WAYLAND)
  file(GLOB GLFW_WAYLAND_SOURCES ${GLFW_SOURCE_DIR}/src/wl_platform.h
                                 ${GLFW_SOURCE_DIR}/src/xkb_unicode.h
                                 ${GLFW_SOURCE_DIR}/src/xkb_unicode.c
                                 ${GLFW_SOURCE_DIR}/src/wl_init.c
                                 ${GLFW_SOURCE_DIR}/src/wl_monitor.c
                                 ${GLFW_SOURCE_DIR}/src/wl_window.c
  )
  target_compile_definitions(glfw PRIVATE _GLFW_WAYLAND)
  target_sources(glfw PRIVATE ${GLFW_WAYLAND_SOURCES})
endif()

if (GLFW_BUILD_X11 OR GLFW_BUILD_WAYLAND)
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_sources(glfw PRIVATE ${GLFW_SOURCE_DIR}/src/linux_joystick.h 
                                ${GLFW_SOURCE_DIR}/src/linux_joystick.c)
  endif()
  target_sources(glfw PRIVATE ${GLFW_SOURCE_DIR}/src/posix_poll.h 
                              ${GLFW_SOURCE_DIR}/src/posix_poll.c)
endif()

if (GLFW_BUILD_WAYLAND)
  include(CheckIncludeFiles)
  include(CheckFunctionExists)
  check_function_exists(memfd_create HAVE_MEMFD_CREATE)
  if (HAVE_MEMFD_CREATE)
    target_compile_definitions(glfw PRIVATE HAVE_MEMFD_CREATE)
  endif()

  find_program(WAYLAND_SCANNER_EXECUTABLE NAMES wayland-scanner)
  if (NOT WAYLAND_SCANNER_EXECUTABLE)
    message(FATAL_ERROR "Failed to find wayland-scanner")
  endif()

  set(GENERATED_WAYLAND_DIR "${CMAKE_BINARY_DIR}/wayland_files")
  file(MAKE_DIRECTORY ${GENERATED_WAYLAND_DIR})
  include_directories(${GENERATED_WAYLAND_DIR})

  macro(generate_wayland_protocol protocol_file)
    set(protocol_path "${GLFW_SOURCE_DIR}/deps/wayland/${protocol_file}")

    string(REGEX REPLACE "\\.xml$" "-client-protocol.h" header_file ${protocol_file})
    string(REGEX REPLACE "\\.xml$" "-client-protocol-code.h" code_file ${protocol_file})

    add_custom_command(
      OUTPUT ${GENERATED_WAYLAND_DIR}/${header_file}
      COMMAND "${WAYLAND_SCANNER_EXECUTABLE}" client-header "${protocol_path}" ${GENERATED_WAYLAND_DIR}/${header_file}
      DEPENDS "${protocol_path}"
      VERBATIM)

    add_custom_command(
      OUTPUT ${GENERATED_WAYLAND_DIR}/${code_file}
      COMMAND "${WAYLAND_SCANNER_EXECUTABLE}" private-code "${protocol_path}" ${GENERATED_WAYLAND_DIR}/${code_file}
      DEPENDS "${protocol_path}"
      VERBATIM)

    target_sources(glfw PRIVATE ${GENERATED_WAYLAND_DIR}/${header_file} ${GENERATED_WAYLAND_DIR}/${code_file})
  endmacro()

  generate_wayland_protocol("wayland.xml")
  generate_wayland_protocol("viewporter.xml")
  generate_wayland_protocol("xdg-shell.xml")
  generate_wayland_protocol("idle-inhibit-unstable-v1.xml")
  generate_wayland_protocol("pointer-constraints-unstable-v1.xml")
  generate_wayland_protocol("relative-pointer-unstable-v1.xml")
  generate_wayland_protocol("fractional-scale-v1.xml")
  generate_wayland_protocol("xdg-activation-v1.xml")
  generate_wayland_protocol("xdg-decoration-unstable-v1.xml")
endif()

target_link_libraries(glfw PRIVATE Threads::Threads)

# Workaround for CMake not knowing about .m files before version 3.16
if (CMAKE_VERSION VERSION_LESS "3.16" AND APPLE)
  set_source_files_properties(${GLFW_SOURCE_DIR}/src/cocoa_joystick.m
                              ${GLFW_SOURCE_DIR}/src/cocoa_init.m
                              ${GLFW_SOURCE_DIR}/src/cocoa_monitor.m
                              ${GLFW_SOURCE_DIR}/src/cocoa_window.m
                              ${GLFW_SOURCE_DIR}/src/nsgl_context.m 
                              PROPERTIES LANGUAGE C
  )
endif()

if (GLFW_BUILD_WIN32)
  list(APPEND glfw_PKG_LIBS "-lgdi32")
endif()

if (GLFW_BUILD_COCOA)
  target_link_libraries(glfw PRIVATE "-framework Cocoa"
                                      "-framework IOKit"
                                      "-framework CoreFoundation"
                                      "-framework QuartzCore")

  set(glfw_PKG_DEPS "")
  set(glfw_PKG_LIBS "-framework Cocoa -framework IOKit -framework CoreFoundation -framework QuartzCore")
endif()

if (GLFW_BUILD_WAYLAND)
  include(FindPkgConfig)

  pkg_check_modules(Wayland REQUIRED
    wayland-client>=0.2.7
    wayland-cursor>=0.2.7
    wayland-egl>=0.2.7
    xkbcommon>=0.5.0)

  target_include_directories(glfw PRIVATE ${Wayland_INCLUDE_DIRS})

  if (NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_package(EpollShim)
    if (EPOLLSHIM_FOUND)
      target_include_directories(glfw PRIVATE ${EPOLLSHIM_INCLUDE_DIRS})
      target_link_libraries(glfw PRIVATE ${EPOLLSHIM_LIBRARIES})
    endif()
  endif()
endif()

if (GLFW_BUILD_X11)
  find_package(X11 REQUIRED)
  target_include_directories(glfw PRIVATE "${X11_X11_INCLUDE_PATH}")

  # Check for XRandR (modern resolution switching and gamma control)
  if (NOT X11_Xrandr_INCLUDE_PATH)
    message(FATAL_ERROR "RandR headers not found; install libxrandr development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xrandr_INCLUDE_PATH}")

  # Check for Xinerama (legacy multi-monitor support)
  if (NOT X11_Xinerama_INCLUDE_PATH)
    message(FATAL_ERROR "Xinerama headers not found; install libxinerama development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xinerama_INCLUDE_PATH}")

  # Check for Xkb (X keyboard extension)
  if (NOT X11_Xkb_INCLUDE_PATH)
    message(FATAL_ERROR "XKB headers not found; install X11 development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xkb_INCLUDE_PATH}")

  # Check for Xcursor (cursor creation from RGBA images)
  if (NOT X11_Xcursor_INCLUDE_PATH)
    message(FATAL_ERROR "Xcursor headers not found; install libxcursor development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xcursor_INCLUDE_PATH}")

  # Check for XInput (modern HID input)
  if (NOT X11_Xi_INCLUDE_PATH)
    message(FATAL_ERROR "XInput headers not found; install libxi development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xi_INCLUDE_PATH}")

  # Check for X Shape (custom window input shape)
  if (NOT X11_Xshape_INCLUDE_PATH)
    message(FATAL_ERROR "X Shape headers not found; install libxext development package")
  endif()
  target_include_directories(glfw PRIVATE "${X11_Xshape_INCLUDE_PATH}")
endif()

if (UNIX AND NOT APPLE)
  find_library(RT_LIBRARY rt)
  mark_as_advanced(RT_LIBRARY)
  if (RT_LIBRARY)
    target_link_libraries(glfw PRIVATE "${RT_LIBRARY}")
    list(APPEND glfw_PKG_LIBS "-lrt")
  endif()

  find_library(MATH_LIBRARY m)
  mark_as_advanced(MATH_LIBRARY)
  if (MATH_LIBRARY)
    target_link_libraries(glfw PRIVATE "${MATH_LIBRARY}")
    list(APPEND glfw_PKG_LIBS "-lm")
  endif()

  if (CMAKE_DL_LIBS)
    target_link_libraries(glfw PRIVATE "${CMAKE_DL_LIBS}")
    list(APPEND glfw_PKG_LIBS "-l${CMAKE_DL_LIBS}")
  endif()
endif()

if (WIN32)
  if (GLFW_USE_HYBRID_HPG)
    target_compile_definitions(glfw PRIVATE _GLFW_USE_HYBRID_HPG)
  endif()
endif()

# Enable a reasonable set of warnings
# NOTE: The order matters here, Clang-CL matches both MSVC and Clang
if (MSVC)
  target_compile_options(glfw PRIVATE "/W3")
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
        CMAKE_C_COMPILER_ID STREQUAL "Clang" OR
        CMAKE_C_COMPILER_ID STREQUAL "AppleClang")

  target_compile_options(glfw PRIVATE "-Wall")
endif()

if (GLFW_BUILD_WIN32)
  target_compile_definitions(glfw PRIVATE UNICODE _UNICODE)
endif()

# HACK: When building on MinGW, WINVER and UNICODE need to be defined before
# the inclusion of stddef.h (by glfw3.h), which is itself included before
# win32_platform.h.  We define them here until a saner solution can be found
# NOTE: MinGW-w64 and Visual C++ do /not/ need this hack.
if (MINGW)
    target_compile_definitions(glfw PRIVATE WINVER=0x0501)
endif()

# Workaround for legacy MinGW not providing XInput and DirectInput
if (MINGW)
    include(CheckIncludeFile)
    check_include_file(dinput.h DINPUT_H_FOUND)
    check_include_file(xinput.h XINPUT_H_FOUND)
    if (NOT DINPUT_H_FOUND OR NOT XINPUT_H_FOUND)
        target_include_directories(glfw PRIVATE "${GLFW_SOURCE_DIR}/deps/mingw")
    endif()
endif()

# Workaround for the MS CRT deprecating parts of the standard library
if (MSVC OR CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
    target_compile_definitions(glfw PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# Workaround for -std=c99 on Linux disabling _DEFAULT_SOURCE (POSIX 2008 and more)
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_compile_definitions(glfw PRIVATE _DEFAULT_SOURCE)
endif() 

foreach(arg ${glfw_PKG_DEPS})
    string(APPEND deps " ${arg}")
endforeach()
foreach(arg ${glfw_PKG_LIBS})
    string(APPEND libs " ${arg}")
endforeach()

set(GLFW_PKG_CONFIG_REQUIRES_PRIVATE "${deps}" CACHE INTERNAL
    "GLFW pkg-config Requires.private")
set(GLFW_PKG_CONFIG_LIBS_PRIVATE "${libs}" CACHE INTERNAL
    "GLFW pkg-config Libs.private")

# Glad library
file(GLOB GLAD_SOURCES ${GLAD_SOURCE_DIR}/src/glad.c)
add_library(glad ${GLAD_SOURCES})

#--------------------------------------------------------------------

set(CGAL_GLFW_FOUND TRUE)
set_property(GLOBAL PROPERTY CGAL_GLFW_FOUND TRUE)

function(CGAL_setup_CGAL_GLFW_dependencies target)
  target_link_libraries(${target} INTERFACE CGAL::CGAL)
  target_link_libraries(${target} INTERFACE glfw)
  target_link_libraries(${target} INTERFACE glad)
 
  # Remove -Wdeprecated-copy for GCC >= 9
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
    target_compile_options(${target} INTERFACE "-Wno-deprecated-copy" "-Wno-cast-function-type")
  endif()
endfunction()
