// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>

#ifndef CGAL_MSVC_STANDARD_HEADER_FIXES_H
#define CGAL_MSVC_STANDARD_HEADER_FIXES_H

#pragma warning(once: 4291)
#pragma warning(once:4503)

#include <cmath>
namespace std {
	using ::sqrt;
}

#include <cstddef>
namespace std{
  using ::size_t;
  using ::ptrdiff_t;
}

#include <ctime>
namespace std{
  using ::time_t;
}

#endif // CGAL_MSVC_STANDARD_HEADER_FIXES_H
