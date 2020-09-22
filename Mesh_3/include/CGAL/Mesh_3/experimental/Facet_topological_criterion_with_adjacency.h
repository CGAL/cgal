// Copyright (c) 2012-2016 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H
#define CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/mesh_standard_criteria.h>
#include <CGAL/number_utils.h>
#include <CGAL/array.h>

namespace CGAL {

namespace Mesh_3 {

template <typename Tr, typename MeshDomain, typename Visitor_>
class Facet_topological_criterion_with_adjacency :
public Mesh_3::Abstract_criterion<Tr, Visitor_>
{
private:
  typedef typename Tr::Facet Facet;

  typedef Mesh_3::Abstract_criterion<Tr,Visitor_> Base;
  typedef typename Base::Quality Quality;
  typedef typename Base::Is_bad  Is_bad;

  typedef Facet_topological_criterion_with_adjacency<Tr,MeshDomain, Visitor_> Self;

  typedef typename Tr::Geom_traits::FT FT;

  const MeshDomain* domain;

public:
  /// Constructor
  Facet_topological_criterion_with_adjacency(const MeshDomain* domain)
    : domain(domain)
  {}

  /// Destructor
  virtual ~Facet_topological_criterion_with_adjacency() {}

protected:
  virtual void do_accept(Visitor_& v) const
  {
    v.visit(*this);
  }

  virtual Self* do_clone() const
  {
    // Call copy ctor on this
    return new Self(*this);
  }

  virtual Is_bad do_is_bad (const Tr& /*tr*/, const Facet& f) const
  {
    typedef typename Tr::Vertex_handle  Vertex_handle;
    typedef typename Tr::Cell_handle    Cell_handle;

    const Cell_handle& ch = f.first;
    const int& i = f.second;

    typedef typename MeshDomain::Surface_patch_index Patch_index;

    const Patch_index patch_index = ch->surface_patch_index(i);

    typedef std::vector<Patch_index> Index_set;

    int nb_vertices_on_curves = 0;

    for(int k = 0; k < 3; ++k) {
      const Vertex_handle v = ch->vertex((i+k+1)&3);
      switch(v->in_dimension()) {
      case 0:
        {
          ++nb_vertices_on_curves; // corners are on curves
          const typename MeshDomain::Corner_index corner_id =
            domain->corner_index(v->index());

          Index_set set;
          domain->get_corner_incidences(corner_id, std::back_inserter(set));
          if(std::find(set.begin(), set.end(), patch_index) == set.end()) {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
            std::cerr << "Bad facet "
                         "(Facet_topological_criterion_with_adjacency: corner #"
                      << corner_id << ", point " << v->point()
                      << ", is not incident to patch #" << patch_index << ")"
                      << std::endl;
#endif
            return Is_bad(Quality(1)); // bad!
          }
        }
        break;
      case 1:
        {
          ++nb_vertices_on_curves;
          const typename MeshDomain::Curve_index curve_id =
            domain->curve_index(v->index());
          Index_set set;
          domain->get_incidences(curve_id, std::back_inserter(set));
          if(std::find(set.begin(), set.end(), patch_index) == set.end()) {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
            std::cerr << "Bad facet "
                         "(Facet_topological_criterion_with_adjacency: curve #"
                      << curve_id << ", at point " << v->point()
                      << ", is not incident to patch #" << patch_index << ")"
                      << std::endl;
#endif
            return Is_bad(Quality(1)); // bad!
          }
        }
        break;
      case 2:
        if(domain->surface_patch_index(v->index()) != patch_index) {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
          std::cerr << "Bad facet (Facet_topological_criterion_with_adjacency: "
                       "vertex at point "
                    << v->point() << " is not on patch #" << patch_index << ")"
                    << std::endl;
#endif
          return Is_bad(Quality(1)); // bad!
        }
        break;
      default:
        return Is_bad(Quality(1));
        break;
      }
    }
    if(nb_vertices_on_curves == 3) {
#ifdef CGAL_MESH_3_DEBUG_FACET_CRITERIA
      std::cerr << "Bad facet (Facet_topological_criterion_with_adjacency: "
                   "three points on a curve)\n";
#endif
      return Is_bad(Quality(1)); // bad!
      // All vertices are on curves. That means that the facet could be on
      // several different patches. Let's disallow that.
    }
    return Is_bad();
  }
}; // end class Facet_topological_criterion_with_adjacency

} // end namespace Mesh_3

} // end namespace CGAL

#endif // CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H
