// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SURFACE_MESH_DEFAULT_CRITERIA_3_H
#define CGAL_SURFACE_MESH_DEFAULT_CRITERIA_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/Surface_mesher/Standard_criteria.h>
#include <iostream>

namespace CGAL {

template <class Tr>
class Surface_mesh_default_criteria_3
{
  typedef Surface_mesher::Refine_criterion<Tr> Criterion;
  typedef Surface_mesher::Standard_criteria<Criterion> Criteria;

public:
  typedef Tr Triangulation;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename Criteria::Quality Quality;
  typedef typename Tr::Facet Facet;

  Surface_mesh_default_criteria_3(const FT angle_bound,
                                  const FT radius_bound,
                                  const FT distance_bound)
    : curvature_size_criterion(distance_bound),
      uniform_size_criterion(radius_bound),
      aspect_ratio_criterion(angle_bound)

  {
    criterion_vector.reserve(4);

    criterion_vector.push_back (&aspect_ratio_criterion);
    criterion_vector.push_back (&uniform_size_criterion);
    criterion_vector.push_back (&curvature_size_criterion);

    criteria.set_criteria(criterion_vector);
  }

  bool is_bad (const Facet& f, Quality& q) const
  {
    return criteria.is_bad(f, q);
  }
private:
  Surface_mesher::Curvature_size_criterion<Tr> curvature_size_criterion;
  // bound on Hausdorff distance does not play any role if bigger than
  // the square of the Uniform_size_criterion

  Surface_mesher::Uniform_size_criterion<Tr> uniform_size_criterion;
  // bound on radii of surface Delaunay balls

  Surface_mesher::Aspect_ratio_criterion<Tr> aspect_ratio_criterion;
  // lower bound on minimum angle in degrees

  std::vector<Criterion*> criterion_vector;
  Criteria criteria;
}; // end class Surface_mesh_default_criteria_3

template <typename Tr>
std::ostream&
operator<<(std::ostream& os,
           const typename Surface_mesh_default_criteria_3<Tr>::Quality& /*q*/)
{
  return os << "q";
}

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_DEFAULT_CRITERIA_3_H
