// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MARCHING_TETRAHEDRA_H
#define CGAL_MARCHING_TETRAHEDRA_H

#include <CGAL/license/Skin_surface_3.h>


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>

namespace CGAL { 

// If TriangulationDataStructure_3 only gets a Cell_handle range
// it is not possible to derive the Vertex_handle type.
template <class Vertex_iterator,
	  class Cell_iterator,
	  class HalfedgeDS_3,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3 >
class Marching_tetrahedra_builder : public Modifier_base<HalfedgeDS_3> {
public:
  typedef HalfedgeDS_3                          HDS;
  typedef MarchingTetrahedraTraits_3            Traits;
  typedef MarchingTetrahedraObserver_3          Observer;

private:
  typedef Vertex_iterator                               T_Vertex_iterator;
  typedef Cell_iterator                                 T_Cell_iterator;
  typedef typename T_Cell_iterator::value_type          T_Cell;
  typedef typename T_Cell::Vertex_handle                T_Vertex_handle;
  typedef typename T_Vertex_iterator::value_type        T_Vertex;
  typedef typename T_Vertex_iterator::value_type *      T_Vertex_pointer;
  typedef typename T_Cell_iterator::value_type *        T_Cell_pointer;

  typedef typename HDS::Face_handle                     HDS_face_handle;
  typedef typename HDS::Vertex_handle                   HDS_vertex_handle;

  typedef std::map<T_Vertex_handle,bool>                T_vertex_map;
  typedef typename T_vertex_map::iterator               T_vertex_map_it;

  // First vertex lies inside the surface, the second vertex outside
  typedef std::pair<T_Vertex_handle,T_Vertex_handle>    T_edge;
  typedef std::map<T_edge,int>                          T_edge_map;
  typedef typename T_edge_map::iterator                 T_edge_map_it;

  typedef Polyhedron_incremental_builder_3<HDS>         Polyh_incr_builder;
public:
  
  Marching_tetrahedra_builder(
    T_Vertex_iterator vertices_begin,
    T_Vertex_iterator vertices_end,
    T_Cell_iterator cells_begin,
    T_Cell_iterator cells_end,
    const Traits &traits, 
    Observer &observer)
    : vertices_begin(vertices_begin),
      vertices_end(vertices_end),
      cells_begin(cells_begin),
      cells_end(cells_end),
      traits(traits), observer(observer), nVertices(0) {
  }

  void operator()( HDS& hds) {
    // sortedV is an array of vertex indices.
    // The first nIn vertices lie inside the surface, the last 4-nIn outside.
    int sortedV[4], cellV[4], nIn;
    T_edge edge;

    Polyh_incr_builder B( hds, true);
    B.begin_surface(0,0,0);

    for (T_Cell_iterator cit = cells_begin; cit != cells_end; ++cit) {
      // Compute signs on vertices and sort them:
      nIn = 0;
      for (int i=0; i<4; i++) {
        if (is_inside(cit,i)) {
          sortedV[nIn] = i; nIn++;
        } else {
          sortedV[3-i+nIn] = i;
        }
      }

      // Process edges whose vertices lie on different sides of the surface
      int edgeN=0;
      for (int i=0; i<nIn; i++) {
        for (int j=nIn; j<4; j++) {
          cellV[edgeN++] = process_edge(B, cit,sortedV[i],sortedV[j]);
        }
      }

      // Construct triangles:
      CGAL_assertion(!B.error());
      if (nIn==1) {
	process_cell(B, cellV, (sortedV[0]%2)==1, cit);
	CGAL_assertion(!B.error());
      } else if (nIn==2) {
	bool change_orientation =
	  (((sortedV[0] == 0) && (sortedV[1]==1)) || 
	    ((sortedV[0] == 0) && (sortedV[1]==3)) || 
	    ((sortedV[0] == 1) && (sortedV[1]==2)) || 
	    ((sortedV[0] == 2) && (sortedV[1]==3)));
	
	process_cell(B, cellV, !change_orientation, cit);
 	process_cell(B, cellV+1, change_orientation, cit);
	CGAL_assertion(!B.error());
      } else if (nIn==3) {
	process_cell(B, cellV, (sortedV[3]%2) == 1, cit);
	CGAL_assertion(!B.error());
      }
      
    }

    CGAL_assertion(!B.error());
    B.end_surface();
  }


  bool is_inside(T_Cell_iterator ch, int i) {
    CGAL_assertion(&*ch != NULL);
    //return (traits.sign(ch,i) == POSITIVE);
    T_vertex_map_it it = triang_vertex_signs.find((ch->vertex(i)));
    
    if (it == triang_vertex_signs.end()) {
      bool side = (traits.sign(ch,i) == POSITIVE);
      CGAL_assertion(triang_vertex_signs.find((ch->vertex(i))) ==
		     triang_vertex_signs.end());
      CGAL_assertion(&*ch != NULL);
      triang_vertex_signs[(ch->vertex(i))] = side;
      CGAL_assertion(triang_vertex_signs.find((ch->vertex(i))) !=
		     triang_vertex_signs.end());
      CGAL_assertion(triang_vertex_signs[(ch->vertex(i))] == side);
      return side;
    } else {
      return it->second;
    }
  }
  int process_edge(Polyh_incr_builder &B,
    T_Cell_iterator ch, int i, int j) {
    CGAL_assertion(is_inside(ch, i));
    CGAL_assertion(!is_inside(ch, j));
    
    T_edge edge = T_edge(ch->vertex(i),ch->vertex(j));
    T_edge_map_it edge_it = polyh_vert.find(edge);

    if (edge_it == polyh_vert.end()) {
      HDS_vertex_handle vh = B.add_vertex(traits.intersection(ch, i, j));
      polyh_vert[edge] = nVertices;
      nVertices ++;
      observer.after_vertex_insertion(ch, i, j, vh);
      return nVertices-1;
    } else {
      return edge_it->second;
    }
  }
  
  // Orientation is right
  void process_cell(
    Polyh_incr_builder &B,
    int *vs,
    bool change_orientation,
    T_Cell_iterator ch)
  {
    CGAL_assertion((vs[0]!=vs[1]) && (vs[0]!=vs[2]) && (vs[1]!=vs[2]));
    HDS_face_handle f = B.begin_facet();
    if (change_orientation) {
      B.add_vertex_to_facet( vs[0] );
      B.add_vertex_to_facet( vs[2] );
      B.add_vertex_to_facet( vs[1] );
    } else {
      B.add_vertex_to_facet( vs[0] );
      B.add_vertex_to_facet( vs[1] );
      B.add_vertex_to_facet( vs[2] );
    }
    B.end_facet();
    observer.after_facet_insertion(ch, f);
  }

private:
  T_Vertex_iterator vertices_begin, vertices_end;
  T_Cell_iterator cells_begin, cells_end;
  const Traits &traits;
  Observer &observer;
  T_edge_map polyh_vert;
  T_vertex_map triang_vertex_signs;
  int nVertices;
};


template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3>
void marching_tetrahedra_3(
  const Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits) {
  
  typedef Marching_tetrahedra_observer_default_3
    <typename Triangulation_3::Finite_vertices_iterator,
    typename Triangulation_3::Finite_cells_iterator, 
    Polyhedron_3>
                                                                Observer; 

  marching_tetrahedra_3(triangulation.finite_vertices_begin(),
			triangulation.finite_vertices_end(),
			triangulation.finite_cells_begin(),
			triangulation.finite_cells_end(),
			polyhedron,
			traits,
			Observer());
}

template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3>
void marching_tetrahedra_3(
  const Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits,
  MarchingTetrahedraObserver_3 &observer) {

  marching_tetrahedra_3(triangulation.finite_vertices_begin(),
			triangulation.finite_vertices_end(),
			triangulation.finite_cells_begin(),
			triangulation.finite_cells_end(),
			polyhedron,
			traits,
			observer);
}

template <class Vertex_iterator,
	  class Cell_iterator,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3>
void marching_tetrahedra_3(
  Vertex_iterator finite_vertices_begin,
  Vertex_iterator finite_vertices_end,
  Cell_iterator finite_cells_begin,
  Cell_iterator finite_cells_end,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits,
  MarchingTetrahedraObserver_3 &observer) {
  
  typedef typename Polyhedron_3::HalfedgeDS                   HDS;
  typedef Marching_tetrahedra_builder<
    Vertex_iterator, Cell_iterator, HDS,
    MarchingTetrahedraTraits_3, MarchingTetrahedraObserver_3> Builder;
  
  Builder builder(finite_vertices_begin, finite_vertices_end,
		  finite_cells_begin, finite_cells_end,
		  traits, observer);
  polyhedron.delegate(builder);
}

} //namespace CGAL 

#endif // CGAL_MARCHING_TETRAHEDRA_H
