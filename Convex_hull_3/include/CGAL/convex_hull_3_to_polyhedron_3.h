// Copyright (c) 2011  GeometryFactory (France).
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
// 
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H
#define CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H

#include <CGAL/license/Convex_hull_3.h>


#define CGAL_DEPRECATED_HEADER "<CGAL/convex_hull_3_to_polyhedron_3.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/convex_hull_3_to_face_graph.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

namespace CGAL {

template <class HDS,class Triangulation>
class Convex_hull_modifier_from_triangulation_3 : public CGAL::Modifier_base<HDS> {
  typedef std::map<typename Triangulation::Vertex_handle,unsigned> Vertex_map;
  
  const Triangulation& t;
  template <class Builder>
  static unsigned get_vertex_index( Vertex_map& vertex_map,
                                    typename Triangulation::Vertex_handle vh,
                                    Builder& builder,
                                    unsigned& vindex)
  {
    std::pair<typename Vertex_map::iterator,bool>
      res=vertex_map.insert(std::make_pair(vh,vindex));
    if (res.second){
      builder.add_vertex(vh->point());
      ++vindex;
    }
    return res.first->second;
  }
  
public:
  Convex_hull_modifier_from_triangulation_3(const Triangulation& t_):t(t_) 
  {
    CGAL_assertion(t.dimension()==3);
  }
  void operator()( HDS& hds) {
    // Postcondition: `hds' is a valid polyhedral surface.

    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    std::vector<typename Triangulation::Cell_handle>  ch_facets;
    Vertex_map vertex_map;
    t.incident_cells(t.infinite_vertex(),std::back_inserter(ch_facets));
    std::size_t nb_facets=ch_facets.size();
    //start the surface
    B.begin_surface( nb_facets, nb_facets);
    unsigned vindex=0;
    for (typename std::vector<typename Triangulation::Cell_handle>::const_iterator it=
          ch_facets.begin();it!=ch_facets.end();++it)
    {
      unsigned ifv_index= (*it)->index(t.infinite_vertex());
      bool is_even=ifv_index%2==0;
      unsigned i0=get_vertex_index(vertex_map,(*it)->vertex((ifv_index + (is_even?3:1) )%4),B,vindex);
      unsigned i1=get_vertex_index(vertex_map,(*it)->vertex((ifv_index + 2            )%4),B,vindex);
      unsigned i2=get_vertex_index(vertex_map,(*it)->vertex((ifv_index + (is_even?1:3) )%4),B,vindex);
      B.begin_facet();
      B.add_vertex_to_facet( i0 );
      B.add_vertex_to_facet( i1 );
      B.add_vertex_to_facet( i2 );
      B.end_facet();      
    }
    B.end_surface();
  }
};

template<class Triangulation_3,class Polyhedron_3>
CGAL_DEPRECATED void convex_hull_3_to_polyhedron_3(const Triangulation_3& T,Polyhedron_3& P){
  P.clear();
  Convex_hull_modifier_from_triangulation_3<typename Polyhedron_3::HalfedgeDS,Triangulation_3> modifier(T);
  P.delegate(modifier);
}

} //namespace CGAL

#endif //CGAL_CONVEX_HULL_3_TO_POLYHEDRON_3_H
