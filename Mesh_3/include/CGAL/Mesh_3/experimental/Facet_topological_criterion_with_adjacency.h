// Copyright (c) 2012-2016 GeometryFactory Sarl (France).
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
  typedef typename Base::Badness Badness;
  
  typedef Facet_topological_criterion_with_adjacency<Tr,MeshDomain, Visitor_> Self;

  typedef typename Tr::Geom_traits::FT FT;

  const MeshDomain* domain;
  
public:
  /// Constructor
  Facet_topological_criterion_with_adjacency(const MeshDomain* domain)
    : domain(domain)
  {};
  /// Destructor
  virtual ~Facet_topological_criterion_with_adjacency() {};
  
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
  
  virtual Badness do_is_bad (const Facet& f) const
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
          if(std::find(set.begin(), set.end(), patch_index) == set.end())
            return Badness(Quality(1)); // bad!
        }
        break;
      case 1: 
        {
          ++nb_vertices_on_curves;
          const typename MeshDomain::Curve_index curve_id =
            domain->curve_index(v->index());
          Index_set set;
          domain->get_incidences(curve_id, std::back_inserter(set));
          if(std::find(set.begin(), set.end(), patch_index) == set.end())
            return Badness(Quality(1)); // bad!
        }
        break;
      case 2:
        if(domain->surface_patch_index(v->index()) != patch_index)
          return Badness(Quality(1)); // bad!
        break;
      default:
        return Badness(Quality(1));
        break;
      }
    }
    if(nb_vertices_on_curves == 3) {
      return Badness(Quality(1)); // bad!
      // All vertices are on curves. That means that the facet could be on
      // several different patches. Let's disallow that.
    }
    return Badness();
  }
}; // end class Facet_topological_criterion_with_adjacency

} // end namespace Mesh_3

} // end namespace CGAL

#endif // CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H
