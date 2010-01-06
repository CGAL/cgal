// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Mesh_facet_criteria_3 class.
//******************************************************************************

#ifndef CGAL_MESH_FACET_CRITERIA_3_H
#define CGAL_MESH_FACET_CRITERIA_3_H

#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>


namespace CGAL {

template<typename Tr, typename Visitor_ = Mesh_3::Facet_criterion_visitor<Tr> >
class Mesh_facet_criteria_3
{
  typedef Visitor_ Visitor;
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;

  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_facet_criteria_3<Tr> Self;

public:
  typedef typename Visitor::Facet_quality Facet_quality;
  typedef typename Visitor::Facet_badness Facet_badness;

  /**
   * @brief Constructor
   */
  Mesh_facet_criteria_3(const FT& angle_bound,
                        const FT& radius_bound,
                        const FT& distance_bound)
  {
    typedef Mesh_3::Aspect_ratio_criterion<Tr,Visitor> Aspect_criterion;
    typedef Mesh_3::Uniform_size_criterion<Tr,Visitor> Uniform_size_criterion;
    typedef Mesh_3::Curvature_size_criterion<Tr,Visitor> Curvature_criterion;
    typedef Mesh_3::Facet_on_surface_criterion<Tr,Visitor> On_surface_criterion;

    if ( 0 != angle_bound )
      criteria_.add(new Aspect_criterion(angle_bound));

    if ( 0 != radius_bound )
      criteria_.add(new Uniform_size_criterion(radius_bound));

    if ( 0 != distance_bound )
      criteria_.add(new Curvature_criterion(distance_bound));

    criteria_.add(new On_surface_criterion());
  }

  /// Destructor
  ~Mesh_facet_criteria_3() { }

   /**
   * @brief returns the badness of facet \c facet
   * @param facet the facet
   * @return the badness of \c facet
   */
  Facet_badness operator()(const Facet& facet) const
  {
    return criteria_(facet);
  }

private:
  Criteria criteria_;
};  // end class Mesh_facet_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_FACET_CRITERIA_3_H
