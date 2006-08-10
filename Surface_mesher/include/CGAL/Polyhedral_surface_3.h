// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: 
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POLYHEDRAL_SURFACE_3_H
#define CGAL_POLYHEDRAL_SURFACE_3_H

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Data_structure_using_octree_3.h>
#include <CGAL/Surface_mesher/Polyhedral_oracle.h>
#include <iostream>
#include <vector>

namespace CGAL {

template <class GT>
class Polyhedral_surface_3
{
public:
  typedef GT Geom_traits;

  class Normalized_geom_traits : public Geom_traits 
  {
  public:
    typedef typename 
    Kernel_traits<typename Geom_traits::Point_3>::Kernel::Point_3 Point_3;
  };

  typedef Data_structure_using_octree_3<Normalized_geom_traits> Subfacets_octree;
  typedef typename GT::Point_3 Point_3;

  typedef Polyhedral_surface_3<GT> Self;

  typedef Surface_mesher::Polyhedral_oracle<Self> Surface_mesher_traits_3;

  typedef typename Subfacets_octree::Bbox Bbox;

  Polyhedral_surface_3(std::istream& input_file)
    : subfacets_octree(), input_points()
  {
    subfacets_octree.input(input_file,
                           std::back_inserter(input_points));
  }

  Bbox bbox() const
  {
    return subfacets_octree.bbox();
  }

public:
  Subfacets_octree subfacets_octree;
  std::vector<Point_3> input_points;
};

} // end namespace CGAL

#endif // CGAL_POLYHEDRAL_SURFACE_3_H
