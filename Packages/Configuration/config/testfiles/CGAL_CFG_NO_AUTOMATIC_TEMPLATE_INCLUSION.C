// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : various

// CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler whether it is able
// locate template implementation files with suffix .C  automatically.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| When template implementation files are not included in the source files,
//| a compiler may attempt to find the unincluded template bodies
//| automatically. For example, suppose that the following conditions are
//| all true.
//|
//| - template entity ABC::f is declared in file xyz.h
//| - an instantiation of ABC::f is required in a compilation
//| - no definition of ABC::f appears in the source code processed by the
//|   compilation
//| 
//| In this case, the compiler may look to see if the source file xyz.n exists,
//| where n is .c, .C, .cpp, .CPP, .cxx, .CXX, or .cc. If this feature is
//| missing, the flag CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION is set.

#include "template.h"

int main()
{
  A<float> a;
  float f = a.square(1.2);
  (void) f;

  return 0;
}

// EOF //
