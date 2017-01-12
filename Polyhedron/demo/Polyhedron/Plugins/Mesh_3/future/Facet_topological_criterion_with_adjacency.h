// Copyright (c) 2012 GeometryFactory Sarl (France)
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
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H
#define CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H

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

    Vertex_handle v[3];
    for(int k = 0; k < 3; ++k) {
      v[k] = ch->vertex((i+k+1)&3);
    }

    if(v[0]->in_dimension() > 1 
       && v[1]->in_dimension() > 1 
       && v[2]->in_dimension() > 1) return Badness();

    typedef typename MeshDomain::Surface_patch_index Patch_index;
    typedef std::vector<Patch_index> Index_set;

    CGAL::cpp0x::array<Index_set,3> sets;

    for(int k = 0; k < 3; ++k) {
      switch(v[k]->in_dimension()) {
      case 0:
        {
          typename MeshDomain::Corner_index corner_id = 
            domain->corner_index(v[k]->index());

          domain->get_corner_incidences(corner_id,
                                        std::inserter(sets[k],
                                                      sets[k].end()));
        }
        break;
      case 1: 
        {
          typename MeshDomain::Curve_segment_index curve_id =
            domain->curve_segment_index(v[k]->index());
          domain->get_incidences(curve_id,
                                 std::inserter(sets[k],
                                               sets[k].end()));
        }
        break;
      case 2:
        sets[k].push_back(domain->surface_patch_index(v[k]->index()));
        break;
      default:
        return Badness(Quality(1));
        break;
      }
    }

    Index_set v0_v1;
    std::set_intersection(sets[0].begin(),
                          sets[0].end(),
                          sets[1].begin(),
                          sets[1].end(),
                          std::inserter(v0_v1, v0_v1.end()));
    Index_set v0_v1_v2;
    std::set_intersection(sets[2].begin(),
                          sets[2].end(),
                          v0_v1.begin(),
                          v0_v1.end(),
                          std::inserter(v0_v1_v2, v0_v1_v2.end()));
    if(v0_v1_v2.empty())
      return Badness(Quality(1)); // bad!
    else
      return Badness();
  }
}; // end class Facet_topological_criterion_with_adjacency

} // end namespace Mesh_3

} // end namespace CGAL

#endif // CGAL_MESH_3_FACET_TOPOLOGICAL_CRITERION_WITH_ADJACENCY_H
