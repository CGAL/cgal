// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/mesh_standard_facet_criteria.h>
#include <CGAL/Mesh_facet_topology.h>
#include <CGAL/Mesh_3/Is_mesh_domain_field_3.h>

namespace CGAL {

/*!
  \ingroup PkgMesh3MeshClasses

  The class `Mesh_facet_criteria_3` is a model of `MeshFacetCriteria_3`.
  It provides a uniform bound for the shape criterion,
  a uniform or variable sizing field
  for the size criterion and/or
  a uniform or variable distance field
  for the approximation error criterion.

  \tparam Tr must be identical to the nested type
  `Triangulation` of the instance used as model of
  `MeshComplex_3InTriangulation_3`.

  \cgalModels `MeshFacetCriteria_3`

  \sa `CGAL::Mesh_facet_topology`
  \sa `MeshCriteria_3`
  \sa `MeshFacetCriteria_3`
  \sa `CGAL::Mesh_criteria_3<Tr>`
  \sa `MeshDomainField_3`
  \sa `CGAL::make_mesh_3()`

*/
template<typename Tr,
         typename Visitor_ = Mesh_3::Facet_criterion_visitor_with_features<Tr> >
class Mesh_facet_criteria_3
{
public:
  typedef Visitor_ Visitor;
  typedef typename Visitor::Facet_quality Facet_quality;
  typedef typename Visitor::Is_facet_bad  Is_facet_bad;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor> Abstract_criterion;
private:
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;

  typedef typename Tr::Facet Facet;
  typedef typename Tr::Geom_traits::FT FT;

  typedef Mesh_facet_criteria_3<Tr> Self;

public:
  typedef CGAL::Tag_true Has_manifold_criterion;

#ifdef DOXYGEN_RUNNING
/// \name Types
/// @{

/*!
Numerical type
@todo:  In the code this typedef is private
*/
typedef Tr::FT FT;

/// @}
#endif


  /// \name Creation
  /// @{

  /*!
    Returns an object to serve as criteria for facets.
    \param angle_bound is the lower bound for the angle in degrees of the
    surface mesh facets.
    \param radius_bound is a uniform upper bound
    for the radius of the surface Delaunay balls.
    \param distance_bound is an upper bound for the center-center distances
    of the surface mesh facets.
    \param topology is the set of topological constraints
    which have to be verified by each surface facet. See
    section \ref Mesh_3DelaunayRefinement for further details.
    Note that if one parameter is set to 0, then its corresponding criteria is ignored.
  */
  template < typename Sizing_field, typename Sizing_field2 >
  Mesh_facet_criteria_3(const FT& angle_bound,
                        const Sizing_field & radius_bound,
                        const Sizing_field2& distance_bound,
                        const Mesh_facet_topology topology =
                          FACET_VERTICES_ON_SURFACE)
  {
    if ( FT(0) != angle_bound )
      init_aspect(angle_bound);

    init_radius(radius_bound,
                Mesh_3::Is_mesh_domain_field_3<Tr, Sizing_field>());

    init_distance(distance_bound,
                  Mesh_3::Is_mesh_domain_field_3<Tr, Sizing_field2>());

    init_topo(topology);
  }
#ifdef DOXYGEN_RUNNING
  /*!
    Returns an object to serve as criteria for facets. The types `SizingField` and
    `DistanceField` must
    be models of the concept `MeshDomainField_3`. The behavior and semantic of the arguments are the same
    as above, except that the radius and distance bound parameters are
    functionals instead of constants.

    @todo This one is only in the doc
  */
  template <class SizingField, class DistanceField>
  Mesh_facet_criteria_3(const FT& angle_bound,
                        const SizingField& radius_bound,
                        const DistanceField& distance_bound,
                        Mesh_facet_topology topology = FACET_VERTICES_ON_SURFACE);
#endif
  /// @}


/// Destructor
  ~Mesh_facet_criteria_3() { }

   /**
   * @brief returns whether the facet `facet` is bad or not.
   * @param tr the triangulation within which `facet` lives
   * @param facet the facet
   */
  Is_facet_bad operator()(const Tr& tr, const Facet& facet) const
  {
    return criteria_(tr, facet);
  }

  void add(Abstract_criterion* criterion)
  {
    criteria_.add(criterion);
  }

  Mesh_facet_topology topology() const {
    return topology_;
  }

private:
  void init_aspect(const FT& angle_bound)
  {
    typedef Mesh_3::Aspect_ratio_criterion<Tr,Visitor> Aspect_criterion;
    criteria_.add(new Aspect_criterion(angle_bound));
  }

  void init_radius(const FT& radius_bound, Tag_false)
  {
    if(FT(0) == radius_bound) return;
    typedef Mesh_3::Uniform_size_criterion<Tr,Visitor> Uniform_size_criterion;
    criteria_.add(new Uniform_size_criterion(radius_bound));
  }

  template <typename Sizing_field>
  void init_radius(const Sizing_field& radius_bound, Tag_true)
  {
    typedef Mesh_3::Variable_size_criterion<Tr,Visitor,Sizing_field> Variable_size_criterion;
    criteria_.add(new Variable_size_criterion(radius_bound));
  }

  void init_distance(const FT& distance_bound, Tag_false)
  {
    if(FT(0) == distance_bound) return;
    typedef Mesh_3::Uniform_curvature_size_criterion<Tr,Visitor> Criterion;
    criteria_.add(new Criterion(distance_bound));
  }

  template <typename Sizing_field>
  void init_distance(const Sizing_field& distance_bound, Tag_true)
  {
    typedef Mesh_3::Variable_curvature_size_criterion<Tr,
                                                      Visitor,
                                                      Sizing_field> Criterion;
    criteria_.add(new Criterion(distance_bound));
  }

  void init_topo(const Mesh_facet_topology topology)
  {
    topology_ = topology;
    switch ( topology % MANIFOLD )
    {
      case FACET_VERTICES_ON_SURFACE:
      {
        typedef Mesh_3::Facet_on_surface_criterion<Tr,Visitor> On_surface_criterion;
        criteria_.add(new On_surface_criterion());
        break;
      }

      case FACET_VERTICES_ON_SAME_SURFACE_PATCH:
      case FACET_VERTICES_ON_SAME_SURFACE_PATCH_WITH_ADJACENCY_CHECK:
        // @TODO: Implement adjacency check !
      {
        typedef Mesh_3::Facet_on_same_surface_criterion<Tr,Visitor> Same_surface_criterion;
        criteria_.add(new Same_surface_criterion());
        break;
      }
    }
  }

private:
  Criteria criteria_;
  Mesh_facet_topology topology_;
};  // end class Mesh_facet_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_FACET_CRITERIA_3_H
