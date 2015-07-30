// This header file is a copy of "boost/config/auto_link.hpp" 
// from boost version 1.44.0
// but slightly modified to accomodate CGAL libraries.

// Before CGAL-4.7-beta1, it has been synchronized with
// libs/config/ version boost-1.58.0-39-g15d56c9, file
// include/boost/config/auto_link.hpp

//  (C) Copyright John Maddock 2003.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         auto_link.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Automatic library inclusion for Borland/Microsoft compilers.
  */

/*************************************************************************

USAGE:
~~~~~~

Before including this header you must define one or more of define the following macros:

CGAL_LIB_NAME:           Required: A string containing the basename of the library,
                          for example boost_regex.
CGAL_LIB_TOOLSET:        Optional: the base name of the toolset.
CGAL_BUILD_SHARED_LIBS:  Optional: when set link to dll rather than static library.
CGAL_LIB_DIAGNOSTIC:     Optional: when set the header will print out the name
                          of the library selected (useful for debugging).
CGAL_AUTO_LINK_NOMANGLE: Specifies that we should link to CGAL_LIB_NAME.lib,
                          rather than a mangled-name version.
CGAL_AUTO_LINK_TAGGED:   Specifies that we link to libraries built with the --layout=tagged option.
                          This is essentially the same as the default name-mangled version, but without
                          the compiler name and version, or the Boost version.  Just the build options.

These macros will be undef'ed at the end of the header, further this header
has no include guards - so be sure to include it only once from your library!

Algorithm:
~~~~~~~~~~

Libraries for Borland and Microsoft compilers are automatically
selected here, the name of the lib is selected according to the following
formula:

CGAL_LIB_PREFIX
   + CGAL_LIB_NAME
   + "_"
   + CGAL_LIB_TOOLSET
   + CGAL_LIB_THREAD_OPT
   + CGAL_LIB_RT_OPT
   "-"
   + CGAL_VERSION

These are defined as:

CGAL_LIB_PREFIX:     "lib" for static libraries otherwise "".

CGAL_LIB_NAME:       The base name of the lib ( for example boost_regex).

CGAL_LIB_TOOLSET:    The compiler toolset name (vc71, vc80 etc).

CGAL_LIB_THREAD_OPT: "-mt" for multithread builds, otherwise nothing.

CGAL_LIB_RT_OPT:     A suffix that indicates the runtime library used,
                      contains one or more of the following letters after
                      a hyphen:

                      s      static runtime (dynamic if not present).
                      g      debug/diagnostic runtime (release if not present).
                      y      Python debug/diagnostic runtime (release if not present).
                      d      debug build (release if not present).
                      p      STLport build.
                      n      STLport build without its IOStreams.

CGAL_VERSION:        Defined in <CGAL/version.h>


***************************************************************************/

#ifndef CGAL_NO_AUTOLINK

#ifdef __cplusplus
#  ifndef BOOST_CONFIG_HPP
#     include <boost/config.hpp>
#  endif
#elif defined(_MSC_VER) && !defined(__MWERKS__) && !defined(__EDG_VERSION__)
//
// C language compatability (no, honestly)
//
#  define BOOST_MSVC _MSC_VER
#  define BOOST_STRINGIZE(X) BOOST_DO_STRINGIZE(X)
#  define BOOST_DO_STRINGIZE(X) #X
#endif
//
// Only include what follows for known and supported compilers:
//
#if defined(BOOST_MSVC) \
    || defined(__BORLANDC__) \
    || (defined(__MWERKS__) && defined(_WIN32) && (__MWERKS__ >= 0x3000)) \
    || (defined(__ICL) && defined(_MSC_EXTENSIONS) && (_MSC_VER >= 1200))

#ifndef CGAL_VERSION
#  include <CGAL/version.h>
#endif

#ifndef CGAL_LIB_NAME
#  error "Macro CGAL_LIB_NAME not set (internal error)"
#endif

//
// error check:
//
#if defined(__MSVC_RUNTIME_CHECKS) && !defined(_DEBUG)
#  pragma message("Using the /RTC option without specifying a debug runtime will lead to linker errors")
#  pragma message("Hint: go to the code generation options and switch to one of the debugging runtimes")
#  error "Incompatible build options"
#endif
//
// select toolset if not defined already:
//
#ifndef CGAL_LIB_TOOLSET
#  if defined(BOOST_MSVC) && (BOOST_MSVC < 1200)
    // Note: no compilers before 1200 are supported
#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1300)

#    ifdef UNDER_CE
       // eVC4:
#      define CGAL_LIB_TOOLSET "evc4"
#    else
       // vc6:
#      define CGAL_LIB_TOOLSET "vc6"
#    endif

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1310)

     // vc7:
#    define CGAL_LIB_TOOLSET "vc7"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1400)

     // vc71:
#    define CGAL_LIB_TOOLSET "vc71"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1500)

     // vc80:
#    define CGAL_LIB_TOOLSET "vc80"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1600)

     // vc90:
#    define CGAL_LIB_TOOLSET "vc90"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1700)

     // vc10:
#    define CGAL_LIB_TOOLSET "vc100"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1800)

     // vc11:
#    define CGAL_LIB_TOOLSET "vc110"

#  elif defined(BOOST_MSVC) && (BOOST_MSVC < 1900)

     // vc12:
#    define CGAL_LIB_TOOLSET "vc120"

# elif defined(BOOST_MSVC)

   // vc14:
#  define CGAL_LIB_TOOLSET "vc140"

#  elif defined(__BORLANDC__)

     // CBuilder 6:
#    define CGAL_LIB_TOOLSET "bcb"

#  elif defined(__ICL)

     // Intel C++, no version number:
#    define CGAL_LIB_TOOLSET "iw"

#  elif defined(__MWERKS__) && (__MWERKS__ <= 0x31FF )

     // Metrowerks CodeWarrior 8.x
#    define CGAL_LIB_TOOLSET "cw8"

#  elif defined(__MWERKS__) && (__MWERKS__ <= 0x32FF )

     // Metrowerks CodeWarrior 9.x
#    define CGAL_LIB_TOOLSET "cw9"

#  endif
#endif // CGAL_LIB_TOOLSET

//
// select thread opt:
//
#if defined(_MT) || defined(__MT__)
#  define CGAL_LIB_THREAD_OPT "-mt"
#else
#  define CGAL_LIB_THREAD_OPT
#endif

#if defined(_MSC_VER) || defined(__MWERKS__)

#  ifdef _DLL

#     if (defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)) && (defined(_STLP_OWN_IOSTREAMS) || defined(__STL_OWN_IOSTREAMS))

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-gydp"
#        elif defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-gdp"
#        elif defined(_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-gydp"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        elif defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gdp"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-p"
#        endif

#     elif defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-gydpn"
#        elif defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-gdpn"
#        elif defined(_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-gydpn"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        elif defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gdpn"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-pn"
#        endif

#     else

#        if defined(_DEBUG) && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-gyd"
#        elif defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gd"
#        else
#            define CGAL_LIB_RT_OPT
#        endif

#     endif

#  else

#     if (defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)) && (defined(_STLP_OWN_IOSTREAMS) || defined(__STL_OWN_IOSTREAMS))

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-sgydp"
#        elif defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-sgdp"
#        elif defined(_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#             define CGAL_LIB_RT_OPT "-sgydp"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        elif defined(_DEBUG)
#             define CGAL_LIB_RT_OPT "-sgdp"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-sp"
#        endif

#     elif defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#            define CGAL_LIB_RT_OPT "-sgydpn"
#        elif defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-sgdpn"
#        elif defined(_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#             define CGAL_LIB_RT_OPT "-sgydpn"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        elif defined(_DEBUG)
#             define CGAL_LIB_RT_OPT "-sgdpn"
#            pragma message("warning: STLport debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-spn"
#        endif

#     else

#        if defined(_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#             define CGAL_LIB_RT_OPT "-sgyd"
#        elif defined(_DEBUG)
#             define CGAL_LIB_RT_OPT "-sgd"
#        else
#            define CGAL_LIB_RT_OPT "-s"
#        endif

#     endif

#  endif

#elif defined(__BORLANDC__)

//
// figure out whether we want the debug builds or not:
//
#if __BORLANDC__ > 0x561
#pragma defineonoption CGAL_BORLAND_DEBUG -v
#endif
//
// sanity check:
//
#if defined(__STL_DEBUG) || defined(_STLP_DEBUG)
#error "Pre-built versions of the CGAL libraries are not provided in STLport-debug form"
#endif

#  ifdef _RTLDLL

#     if defined(CGAL_BORLAND_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#         define CGAL_LIB_RT_OPT "-yd"
#     elif defined(CGAL_BORLAND_DEBUG)
#         define CGAL_LIB_RT_OPT "-d"
#     elif defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#         define CGAL_LIB_RT_OPT -y
#     else
#         define CGAL_LIB_RT_OPT
#     endif

#  else

#     if defined(CGAL_BORLAND_DEBUG)\
               && defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#         define CGAL_LIB_RT_OPT "-syd"
#     elif defined(CGAL_BORLAND_DEBUG)
#         define CGAL_LIB_RT_OPT "-sd"
#     elif defined(CGAL_DEBUG_PYTHON) && defined(CGAL_LINKING_PYTHON)
#         define CGAL_LIB_RT_OPT "-sy"
#     else
#         define CGAL_LIB_RT_OPT "-s"
#     endif

#  endif

#endif

//
// select linkage opt:
//
#if (defined(_DLL) || defined(_RTLDLL)) && defined(CGAL_BUILD_SHARED_LIBS)
#  define CGAL_LIB_PREFIX
#elif defined(CGAL_BUILD_SHARED_LIBS)
#  error "Mixing a dll CGAL library with a static runtime is a really bad idea..."
#else
#  define CGAL_LIB_PREFIX "lib"
#endif

//
// now include the lib:
//
#if defined(CGAL_LIB_NAME) \
      && defined(CGAL_LIB_PREFIX) \
      && defined(CGAL_LIB_TOOLSET) \
      && defined(CGAL_LIB_THREAD_OPT) \
      && defined(CGAL_LIB_RT_OPT) \
      && defined(CGAL_VERSION)

#ifdef CGAL_AUTO_LINK_TAGGED
#  pragma commentcomment(lib, CGAL_LIB_PREFIX BOOST_STRINGIZE(CGAL_LIB_NAME) CGAL_LIB_THREAD_OPT CGAL_LIB_RT_OPT ".lib")
#  ifdef CGAL_LIB_DIAGNOSTIC
#     pragma message ("Linking to lib file: " CGAL_LIB_PREFIX BOOST_STRINGIZE(CGAL_LIB_NAME) "-" CGAL_LIB_TOOLSET CGAL_LIB_THREAD_OPT CGAL_LIB_RT_OPT "-" BOOST_STRINGIZE(CGAL_VERSION) ".lib")
#  endif
#elif defined(CGAL_AUTO_LINK_NOMANGLE)
#  pragma comment(lib, BOOST_STRINGIZE(CGAL_LIB_NAME) ".lib")
#  ifdef CGAL_LIB_DIAGNOSTIC
#     pragma message ("Linking to lib file: " BOOST_STRINGIZE(CGAL_LIB_NAME) ".lib")
#  endif
#else
#  pragma comment(lib, CGAL_LIB_PREFIX BOOST_STRINGIZE(CGAL_LIB_NAME) "-" CGAL_LIB_TOOLSET CGAL_LIB_THREAD_OPT CGAL_LIB_RT_OPT "-" BOOST_STRINGIZE(CGAL_VERSION) ".lib")
#  ifdef CGAL_LIB_DIAGNOSTIC
#     pragma message ("Linking to lib file: " CGAL_LIB_PREFIX BOOST_STRINGIZE(CGAL_LIB_NAME) "-" CGAL_LIB_TOOLSET CGAL_LIB_THREAD_OPT CGAL_LIB_RT_OPT "-" BOOST_STRINGIZE(CGAL_VERSION) ".lib")
#  endif
#endif

#else
#  error "some required macros where not defined (internal logic error)."
#endif


#endif // _MSC_VER || __BORLANDC__

#endif // CGAL_NO_AUTOLINK

//
// finally undef any macros we may have set:
//
#ifdef CGAL_LIB_PREFIX
#  undef CGAL_LIB_PREFIX
#endif
#if defined(CGAL_LIB_NAME)
#  undef CGAL_LIB_NAME
#endif
// Don't undef this one: it can be set by the user and should be the 
// same for all libraries:
//#if defined(CGAL_LIB_TOOLSET)
//#  undef CGAL_LIB_TOOLSET
//#endif
#if defined(CGAL_LIB_THREAD_OPT)
#  undef CGAL_LIB_THREAD_OPT
#endif
#if defined(CGAL_LIB_RT_OPT)
#  undef CGAL_LIB_RT_OPT
#endif
#if defined(CGAL_LIB_LINK_OPT)
#  undef CGAL_LIB_LINK_OPT
#endif
#if defined(CGAL_LIB_DEBUG_OPT)
#  undef CGAL_LIB_DEBUG_OPT
#endif
#if defined(CGAL_DYN_LINK)
#  undef CGAL_DYN_LINK
#endif
#if defined(CGAL_AUTO_LINK_NOMANGLE)
#  undef CGAL_AUTO_LINK_NOMANGLE
#endif
