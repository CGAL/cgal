// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind.
//
// Every use of CGAL requires a license.
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation.
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres_d_configure.h
// package       : Min_sphere_of_spheres_d (1.10)
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Kaspar Fischer
// maintainer    : Kaspar Fischer <fischerk@inf.ethz.ch>
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


#ifndef CGAL_MINIBALL_CONFIGURE
#define CGAL_MINIBALL_CONFIGURE

// Remark: In case you want to fine-tune the code, feel free to change
// the options below.

// Option: Namespace name
//
// Background: By default, all data-structures and routines of the
// package are placed in a namespace called CGAL.  You can
// change this name by altering the following #define.
//
// Default value: CGAL
#define CGAL_MINIBALL_NAMESPACE CGAL

// Option: Assertions
//
// Background: The package contains lots of assertions (i.e. internal
// consistency checks).  For instance, if assertions are enabled then
// the package will complain when you add balls with negative radii.
// If you disable assertions, such tests will not be made (with the
// advantage that the code is slightly faster).  Do *not* disable
// assertions during development.
//
// Default setting: defined
#ifndef CGAL_NO_ASSERTIONS
#define CGAL_MINIBALL_DEBUG
#endif
#ifdef NDEBUG
#undef  CGAL_MINIBALL_DEBUG
#endif

// (You should not have to alter anything below here.)

// If CGAL is not being used, we need to define certain things:
#ifndef CGAL_VERSION
  namespace CGAL_MINIBALL_NAMESPACE {
    struct Tag_true {};
    struct Tag_false {};
  }
  #define CGAL_MINIBALL_NTS 
#else
  #include <CGAL/basic.h>
  #define CGAL_MINIBALL_NTS CGAL_NTS
#endif

// Define some assertion macros used in the code.
#ifdef CGAL_MINIBALL_DEBUG
  #define CGAL_MINIBALL_ASSERT(expr) assert(expr)
  #define CGAL_MINIBALL_DO_DEBUG(expr) expr
#else
  #define CGAL_MINIBALL_ASSERT(expr) ;
  #define CGAL_MINIBALL_DO_DEBUG(expr) ;
#endif

// Currently, we include all code in the header files because most
// compilers don't support exporting templates anyway:
#define CGAL_MINIBALL_NO_TEMPLATE_EXPORT

// Fixes for Metrowerks CodeWarrior 7 compilers on MacOS X:
// Unfortunately, this compiler puts the following symbols into
// the global namespace, so we put them back into std:
// (This fix is not necessary any more for CW 8 or higher.)
#if defined __MWERKS__ && defined __MACH__ && (__MWERKS__ < 0x3000)
#include <cstdlib>
#include <cmath>
namespace std {
  int rand(void) { return ::rand(); }
  double sqrt(const double x) { return ::sqrt(x); }
  double abs(const double x) { return ::abs(x); }
}
#endif

// Fixes for GCC series 2.95.  (This fix is necessary for 2.95.2 at
// least.  But I guess it is needed for any earlier version of the 2.95
// series, too.)  Apparently, GCC doesn't come with a bitset and sstream
// implementation, that's why we include them here.
#if defined(__GNUC__) && __GNUC__==2 && \
  __GNUC_MINOR__==95 && __GNUC_PATCHLEVEL__ <= 2
  #include <CGAL/Min_sphere_of_spheres_d_gcc2-95-2_fix.h>
#else
  #include <bitset>
  #include <sstream>
#endif

#endif // CGAL_MINIBALL_CONFIGURE
