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
// Author(s)     : Sylvain Pion

//| This flag is set if the compiler bugs with std::numeric_limits

#include <cmath>
#include <limits>

bool
is_finite(double d)
{ return (d != std::numeric_limits<double>::infinity()) && (-d != std::numeric_limits<double>::infinity()); }

int main()
{

  double zero = 0;
  double inf = 1. / zero;
  double nan = zero*inf;
  bool b = true;
  b = b && !is_finite(inf);

  (void) nan; // stop warning

  if (!b)
    return -1;
  return 0;
}
