// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin


#ifndef TANGENTIAL_COMPLEX_H
#define TANGENTIAL_COMPLEX_H

#include <CGAL/basic.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>

#include <vector>

namespace CGAL {
  
/// The class Tangential_complex represents a tangential complex
template <
  typename Kernel, 
  int Intrinsic_dimension, 
  typename Tr = Regular_triangulation<Regular_triangulation_euclidean_traits<
                  CGAL::Epick_d<Dimension_tag<Intrinsic_dimension> > > >
>
class Tangential_complex
{
  typedef typename Kernel::Point_d              Point;
  typedef typename Kernel::Vector_d             Vector;
  typedef typename std::vector<Vector>          Tangent_space_base;

public:
  /// Constructor
  Tangential_complex() {}

  /// Destructor
  ~Tangential_complex() {}

private:
  std::vector<Point>              m_points;
  std::vector<Tangent_space_base> m_tangent_spaces;
  std::vector<Tr>                 m_triangulations;

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
