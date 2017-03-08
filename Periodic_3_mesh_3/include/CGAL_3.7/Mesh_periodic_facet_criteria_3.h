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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Mesh_3/include/CGAL/Mesh_facet_criteria_3.h $
// $Id: Mesh_facet_criteria_3.h 52705 2009-10-23 10:27:15Z stayeb $
//
//
// Author(s)     : Mikhail Bogdanov
//
//******************************************************************************
// File Description :
// Mesh_periodic_facet_criteria_3 class.
//******************************************************************************

#ifndef CGAL_MESH_PERIODIC_FACET_CRITERIA_3_H
#define CGAL_MESH_PERIODIC_FACET_CRITERIA_3_H

#include <CGAL/Mesh_3/periodic_mesh_standard_facet_criteria.h>

namespace CGAL {

template<typename Tr, typename Visitor_ = Mesh_3::Facet_criterion_visitor<Tr> >
class Mesh_periodic_facet_criteria_3
{
  typedef Visitor_ Visitor;
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;

  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;
  typedef typename Tr::Iso_cuboid Iso_cuboid;

  typedef Mesh_periodic_facet_criteria_3<Tr> Self;

public:
  typedef typename Visitor::Facet_quality Facet_quality;
  typedef typename Visitor::Facet_badness Facet_badness;

  /**
   * @brief Constructor
   */
  Mesh_periodic_facet_criteria_3(const Iso_cuboid& periodic_domain,
                                 const FT& angle_bound,
                                 const FT& radius_bound,
                                 const FT& distance_bound) :
    helper_tr(periodic_domain)
  {
    initialize(angle_bound, radius_bound, distance_bound);
  }

  /**
   * @brief Constructor
   */
  template<class PMD>
  Mesh_periodic_facet_criteria_3(const PMD& periodic_mesh_domain,
                                 const FT& angle_bound,
                                 const FT& radius_bound,
                                 const FT& distance_bound) :
    helper_tr(periodic_mesh_domain.periodic_cuboid())
  {
    initialize(angle_bound, radius_bound, distance_bound);
  }

  /// Destructor
  ~Mesh_periodic_facet_criteria_3() { }

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

  /// This function can be called by a constructor only
  void initialize(const FT& angle_bound, 
    const FT& radius_bound, 
    const FT& distance_bound)
  {
    typedef Mesh_3::Periodic_mesh_3::Aspect_ratio_criterion<Tr,Visitor> Aspect_criterion;
    typedef Mesh_3::Periodic_mesh_3::Uniform_size_criterion<Tr,Visitor> Uniform_size_criterion;
    typedef Mesh_3::Periodic_mesh_3::Curvature_size_criterion<Tr,Visitor> Curvature_criterion;
    typedef Mesh_3::Periodic_mesh_3::Facet_on_surface_criterion<Tr,Visitor> On_surface_criterion;

    if ( 0 != angle_bound )
      criteria_.add(new Aspect_criterion(helper_tr, angle_bound));

    if ( 0 != radius_bound )
      criteria_.add(new Uniform_size_criterion(helper_tr, radius_bound));

    if ( 0 != distance_bound )
      criteria_.add(new Curvature_criterion(helper_tr, distance_bound));

    criteria_.add(new On_surface_criterion());
  }

  Criteria criteria_;

  Tr helper_tr;
};  // end class Mesh_facet_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_PERIODIC_FACET_CRITERIA_3_H
