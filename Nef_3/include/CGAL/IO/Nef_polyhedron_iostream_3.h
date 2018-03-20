// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYHEDRON_IOSTREAM_3_H
#define CGAL_NEF_POLYHEDRON_IOSTREAM_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_io_parser.h>

namespace CGAL {

template <typename Kernel, typename Items, typename Mark>
std::ostream& operator<<
 (std::ostream& os, Nef_polyhedron_3<Kernel,Items,Mark>& NP)
{
  typedef typename Nef_polyhedron_3<Kernel,Items, Mark>::SNC_structure SNC_structure;
#ifdef CGAL_NEF3_SORT_OUTPUT
  CGAL::SNC_io_parser<SNC_structure> O(os, NP.snc(), true, true);
#else
  CGAL::SNC_io_parser<SNC_structure> O(os, NP.snc(), false, true);
#endif
  O.print();
  return os;
}

template <typename Kernel, typename Items, typename Mark>
std::istream& operator>>
  (std::istream& is, Nef_polyhedron_3<Kernel,Items, Mark>& NP)
{
  typedef typename Nef_polyhedron_3<Kernel,Items,Mark>::SNC_structure SNC_structure;
  CGAL::SNC_io_parser<SNC_structure> I(is, NP.snc());
  I.read();
  NP.pl()->initialize(&NP.snc());
  return is;
}

} //namespace CGAL

#endif //CGAL_NEF_POLYHEDRON_IOSTREAM_3_H
