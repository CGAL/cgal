// This header file is a copy of "boost/config/auto_link.hpp" 
// from boost version 1.33.1
// but slightly modified to accomodate CGAL libraries.
//-------------------------------------------------------------------------------------- 
//
//  (C) Copyright John Maddock 2003.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//-------------------------------------------------------------------------------------- 
// Copyright (C) 2007 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
// Author(s)     : Fernando Cacciola (fernando.cacciola@geometryfactory.com)
/*************************************************************************

USAGE:
~~~~~~

Before including this header you must define one or more of define the following macros:

CGAL_LIB_NAME:           Required: A string containing the basename of the library,
                          for example boost_regex.
CGAL_LIB_TOOLSET:        Optional: the base name of the toolset.
CGAL_LIB_DIAGNOSTIC:     Optional: when set the header will print out the name
                          of the library selected (useful for debugging).
CGAL_AUTO_LINK_NOMANGLE: Specifies that we should link to CGAL_LIB_NAME.lib,
                          rather than a mangled-name version.

ALL these macros will be undef'ed at the end of the header, even though they are defined from the outside.
That means you must always define them before including this.

Further, this header has no include guards because you can reuse it several times with different
macros (arguments) in order to link different libraries.
Be sure to include it only once for each target library!

Algorithm:
~~~~~~~~~~

Libraries for Borland and Microsoft compilers are automatically
selected here, the name of the lib is selected according to the following
formula:

CGAL_LIB_NAME
   + "_"
   + CGAL_LIB_TOOLSET
   + CGAL_LIB_THREAD_OPT
   + CGAL_LIB_RT_OPT

These are defined as:

CGAL_LIB_NAME:       The base name of the lib ( for example boost_regex).

CGAL_LIB_TOOLSET:    The compiler toolset name (vc71, vc80 etc).

CGAL_LIB_THREAD_OPT: "-mt" for multithread builds, otherwise nothing.

CGAL_LIB_RT_OPT:     A suffix that indicates the runtime library used,
                      contains one or more of the following letters after
                      a hiphen:

                      s      static runtime (dynamic if not present).
                      d      debug build (release if not present).
                      g      debug/diagnostic runtime (release if not present).
                      p      STLPort Build.


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
#  define BOOST_STRINGIZEIZE(X) BOOST_DO_STRINGIZE(X)
#  define BOOST_DO_STRINGIZE(X) #X
#endif
//
// Only include what follows for known and supported compilers:
//
#if defined(BOOST_MSVC) \
    || defined(__BORLANDC__) \
    || (defined(__MWERKS__) && defined(_WIN32) && (__MWERKS__ >= 0x3000)) \
    || (defined(__ICL) && defined(_MSC_EXTENSIONS) && (_MSC_VER >= 1200))

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
#if defined(BOOST_MSVC) && (BOOST_MSVC == 1200)

   // vc6:
#  define CGAL_LIB_TOOLSET "vc6"

#elif defined(BOOST_MSVC) && (BOOST_MSVC == 1300)

   // vc7:
#  define CGAL_LIB_TOOLSET "vc7"

#elif defined(BOOST_MSVC) && (BOOST_MSVC == 1310)

   // vc71:
#  define CGAL_LIB_TOOLSET "vc71"

#elif defined(BOOST_MSVC) && (BOOST_MSVC == 1400)

   // vc80:
#  define CGAL_LIB_TOOLSET "vc80"

#elif defined(BOOST_MSVC) && (BOOST_MSVC == 1500)

   // vc90:
#  define CGAL_LIB_TOOLSET "vc90"

#elif defined(BOOST_MSVC) && (BOOST_MSVC >= 1600)

   // vc10:
#  define CGAL_LIB_TOOLSET "vc100"

#elif defined(__BORLANDC__)

   // CBuilder 6:
#  define CGAL_LIB_TOOLSET "bcb"

#elif defined(__ICL)

   // Intel C++, no version number:
#  define CGAL_LIB_TOOLSET "iw"

#elif defined(__MWERKS__) && (__MWERKS__ <= 0x31FF )

   // Metrowerks CodeWarrior 8.x
#  define CGAL_LIB_TOOLSET "cw8"

#elif defined(__MWERKS__) && (__MWERKS__ <= 0x32FF )

   // Metrowerks CodeWarrior 9.x
#  define CGAL_LIB_TOOLSET "cw9"

#endif
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

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-gdp"
#        elif defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gdp"
#            pragma message("warning: STLPort debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-p"
#        endif

#     elif defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-gdpn"
#        elif defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gdpn"
#            pragma message("warning: STLPort debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-pn"
#        endif

#     else

#        if defined(_DEBUG)
#            define CGAL_LIB_RT_OPT "-gd"
#        else
#            define CGAL_LIB_RT_OPT
#        endif

#     endif

#  else

#     if (defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)) && (defined(_STLP_OWN_IOSTREAMS) || defined(__STL_OWN_IOSTREAMS))

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-sgdp"
#        elif defined(_DEBUG)
#             define CGAL_LIB_RT_OPT "-sgdp"
#            pragma message("warning: STLPort debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-sp"
#        endif

#     elif defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)

#        if defined(_DEBUG) && (defined(__STL_DEBUG) || defined(_STLP_DEBUG))
#            define CGAL_LIB_RT_OPT "-sgdpn"
#        elif defined(_DEBUG)
#             define CGAL_LIB_RT_OPT "-sgdpn"
#            pragma message("warning: STLPort debug versions are built with /D_STLP_DEBUG=1")
#            error "Build options aren't compatible with pre-built libraries"
#        else
#            define CGAL_LIB_RT_OPT "-spn"
#        endif

#     else

#        if defined(_DEBUG)
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
#error "Pre-built versions of the CGAL libraries are not provided in STLPort-debug form"
#endif

#  ifdef _RTLDLL

#     ifdef CGAL_BORLAND_DEBUG
#         define CGAL_LIB_RT_OPT "-d"
#     else
#         define CGAL_LIB_RT_OPT
#     endif

#  else

#     ifdef CGAL_BORLAND_DEBUG
#         define CGAL_LIB_RT_OPT "-sd"
#     else
#         define CGAL_LIB_RT_OPT "-s"
#     endif

#  endif

#endif

//
// now include the lib:
//
#if defined(CGAL_LIB_NAME) \
      && defined(CGAL_LIB_TOOLSET) \
      && defined(CGAL_LIB_THREAD_OPT) \
      && defined(CGAL_LIB_RT_OPT) \

#ifndef CGAL_AUTO_LINK_NOMANGLE
#  define CGAL_LIB_FULL_NAME BOOST_STRINGIZE(CGAL_LIB_NAME) "-" CGAL_LIB_TOOLSET CGAL_LIB_THREAD_OPT CGAL_LIB_RT_OPT ".lib"
#else
#  define CGAL_LIB_FULL_NAME BOOST_STRINGIZE(CGAL_LIB_NAME) ".lib"
#endif

#pragma comment(lib, CGAL_LIB_FULL_NAME )
#ifdef CGAL_LIB_DIAGNOSTIC
#   pragma message ("Linking to lib file: " CGAL_LIB_FULL_NAME )
#endif

#else
#  error "some required macros where not defined (internal logic error)."
#endif


#endif // _MSC_VER || __BORLANDC__

#endif // CGAL_NO_AUTOLINK

//
// finally undef any macros we may have set:
//
#if defined(CGAL_LIB_NAME)
#  undef CGAL_LIB_NAME
#endif
#if defined(CGAL_LIB_TOOLSET)
#  undef CGAL_LIB_TOOLSET
#endif
#if defined(CGAL_LIB_THREAD_OPT)
#  undef CGAL_LIB_THREAD_OPT
#endif
#if defined(CGAL_LIB_RT_OPT)
#  undef CGAL_LIB_RT_OPT
#endif
#if defined(CGAL_LIB_LINK_OPT)
#  undef CGAL_LIB_LINK_OPT
#endif
#if defined(CGAL_AUTO_LINK_NOMANGLE)
#  undef CGAL_AUTO_LINK_NOMANGLE
#endif
#if defined(CGAL_LIB_FULL_NAME)
#  undef CGAL_LIB_FULL_NAME
#endif
