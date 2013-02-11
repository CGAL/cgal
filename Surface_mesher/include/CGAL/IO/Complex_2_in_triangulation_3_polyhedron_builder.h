// Copyright (c) 2003-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008       GeometryFactory (France)
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_POLYHEDRON_BUILDER_H
#define CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_POLYHEDRON_BUILDER_H

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

namespace CGAL {

template <class C2T3, class Polyhedron_>
class Complex_2_in_triangulation_3_polyhedron_builder
  : public CGAL::Modifier_base<typename Polyhedron_::HalfedgeDS>
{
public:
  typedef C2T3 C2t3;
  typedef Polyhedron_ Polyhedron;
  typedef typename Polyhedron::HalfedgeDS HDS;

private:
  typedef typename C2T3::Triangulation Tr;

  const C2t3& c2t3;
  typedef CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> Base;

  template <class IBuilder, class Vertex_handle>
  int get_vertex_index(IBuilder& builder, Vertex_handle vh, std::map<Vertex_handle, int>& V, int& inum)
  {
    typedef typename std::map<Vertex_handle, int>::iterator map_iterator;
    std::pair<map_iterator,bool> insert_res = V.insert( std::make_pair(vh,inum) );
    if ( insert_res.second ){
      typename Tr::Point p = static_cast<typename Tr::Point>(vh->point());
      builder.add_vertex(p);
      ++inum;
    }
    return insert_res.first->second;
  }

public:
  Complex_2_in_triangulation_3_polyhedron_builder(const C2t3& c2t3)
    : Base(), c2t3(c2t3)
  {
  }

  void operator()(HDS& hds) {
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Geom_traits::Vector_3 Vector;
    typedef typename Tr::Edge Edge;
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;

    const Tr& tr = c2t3.triangulation();
    CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
    const typename Tr::size_type number_of_facets = c2t3.number_of_facets();
    builder.begin_surface(tr.number_of_vertices(),
			  number_of_facets);
    {
      // Finite vertices coordinates.
      Finite_facets_iterator fit = tr.finite_facets_begin();
      std::set<Facet> oriented_set;
      std::stack<Facet> stack;

      CGAL_assertion_code(typename Tr::size_type nb_facets = 0; )

      while (oriented_set.size() != number_of_facets) {
        while ( fit->first->is_facet_on_surface(fit->second) == false ||
                oriented_set.find(*fit) != oriented_set.end() ||

                oriented_set.find(c2t3.opposite_facet(*fit)) !=
                oriented_set.end() ) {
          ++fit;
        }
        oriented_set.insert(*fit);
        stack.push(*fit);
        while(! stack.empty() ) {
          Facet f = stack.top();
          stack.pop();
          for(int ih = 0 ; ih < 3 ; ++ih) {
            const int i1  = tr.vertex_triple_index(f.second, tr. cw(ih));
            const int i2  = tr.vertex_triple_index(f.second, tr.ccw(ih));
            if( c2t3.face_status(Edge(f.first, i1, i2)) == C2t3::REGULAR ) {
              Facet fn = c2t3.neighbor(f, ih);
              if (oriented_set.find(fn) == oriented_set.end() &&
                  oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
              {
                oriented_set.insert(fn);
                stack.push(fn);
              }
            } // end "if the edge is regular"
          } // end "for each neighbor of f"
        } // end "stack non empty"
      } // end "oriented_set not full"

      // Orients the whole mesh towards outside:
      // - find the facet with max z
      typename std::set<Facet>::const_iterator top_facet = oriented_set.begin();
      for(typename std::set<Facet>::const_iterator fit = oriented_set.begin();
	  fit != oriented_set.end();
	  ++fit)
      {
	double top_z = 
	  (top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0))->point().z()
	 + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1))->point().z()
	 + top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2))->point().z())/3.;
	double z = 
	  (fit->first->vertex(tr.vertex_triple_index(fit->second, 0))->point().z()
	 + fit->first->vertex(tr.vertex_triple_index(fit->second, 1))->point().z()
	 + fit->first->vertex(tr.vertex_triple_index(fit->second, 2))->point().z())/3.;
        if (top_z < z)
          top_facet = fit;
      }
      // - orient the facet with max z towards +Z axis
      Vertex_handle v0 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 0));
      Vertex_handle v1 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 1));
      Vertex_handle v2 = top_facet->first->vertex(tr.vertex_triple_index(top_facet->second, 2));
      Vector normal = cross_product(v1->point()-v0->point(), v2->point()-v1->point());
      const Vector Z(0, 0, 1);
      bool regular_orientation = (Z * normal >= 0);

      //used to set indices of vertices
      std::map<Vertex_handle, int> V;
      int inum = 0;
      
      for(typename std::set<Facet>::const_iterator fit =
	    oriented_set.begin();
	  fit != oriented_set.end();
	  ++fit)
      {
	int indices[3];
	int index = 0;
	for (int i=0; i<3; i++)
	  indices[index++] = get_vertex_index(
            builder, fit->first->vertex(tr.vertex_triple_index(fit->second, i)), V, inum
          );
	builder.begin_facet();
	  builder.add_vertex_to_facet(indices[0]);
	  builder.add_vertex_to_facet(regular_orientation ? indices[1] : indices[2]);
	  builder.add_vertex_to_facet(regular_orientation ? indices[2] : indices[1]);
	builder.end_facet();
	CGAL_assertion_code(++nb_facets);
      }
      CGAL_assertion(nb_facets == number_of_facets);
      // 	for( Finite_facets_iterator fit = tr.finite_facets_begin();
      // 	     fit != tr.finite_facets_end(); ++fit)
      // 	  if ((*fit).first->is_facet_on_surface((*fit).second)==true)
      // 	  {
      // 	    int indices[3];
      // 	    int index = 0;
      // 	    for (int i=0; i<3; i++)
      // 	      std::cerr << ( indices[index++] = V[(*fit).first->vertex(tr.vertex_triple_index(fit->second, i))] ) << ", ";
      // 	    std::cerr << "\n";
      // 	    if( builder.test_facet(indices+0, indices+3) )
      // 	      builder.add_facet(indices+0, indices+3);
      // 	    else
      // 	    {
      // 	      builder.begin_facet();
      // 	      builder.add_vertex_to_facet(indices[2]);
      // 	      builder.add_vertex_to_facet(indices[1]);
      // 	      builder.add_vertex_to_facet(indices[0]);
      // 	      builder.end_facet();
      // 	    }
      // 	    CGAL_assertion_code(++nb_facets);
      // 	  }
    }
    builder.end_surface();
  } // end operator()
}; // end Complex_2_in_triangulation_3_polyhedron_builder

} // end namespace CGAL

#endif  // CGAL_IO_COMPLEX_2_IN_TRIANGULATION_3_POLYHEDRON_BUILDER_H
