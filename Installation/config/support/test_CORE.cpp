// Copyright (c) 2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
//
// Author(s)     : various

// Tests if CORE is available.

#ifdef CGAL_USE_CGAL_CORE

// We haven't compiled it yet, so no way to test...
int main() { return 0; }

#else

#include <CGAL/CORE/CORE.h>
#include <iostream>

int main()
{
  Expr x = 2, y = 3;
  std::cout << x << " * " << y  << " = " << x*y << std::endl;

  // CORE does not have VERSION macros yet (as of december 2005).
  std::cout << "version=unknown" << std::endl;

  return 0;
}

#endif // CGAL_USE_CGAL_CORE
