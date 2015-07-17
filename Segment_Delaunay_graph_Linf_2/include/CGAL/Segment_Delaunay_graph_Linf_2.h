// Copyright (c) 2015  Università della Svizzera italiana.
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
// Author(s)     : Panagiotis Cheilaris, Sandeep Kumar Dey, Evanthia Papadopoulou
//philaris@gmail.com, sandeep.kr.dey@gmail.com, evanthia.papadopoulou@usi.ch

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_H

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/basic.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_vertex_base_2.h>
#include <CGAL/Segment_Delaunay_graph_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/Segment_Delaunay_graph_2/in_place_edge_list.h>
#include <CGAL/Segment_Delaunay_graph_2/edge_list.h>
#include <CGAL/Segment_Delaunay_graph_2/Traits_wrapper_2.h>
#include <CGAL/Segment_Delaunay_graph_2/Constructions_C2.h>
#include <CGAL/Segment_Delaunay_graph_Linf_2/Constructions_C2.h>

#include <CGAL/Iterator_project.h>
#include <CGAL/utility.h>

#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>

#include <boost/iterator/counting_iterator.hpp>

/*
  Conventions:
  ------------
  1. we treat segments as open; the endpoints are separate objects
  2. a segment of length zero is treated as a point
  3. a point is deleted only if it has no segment adjacent to it
  4. when a segment is deleted it's endpoints are not deleted
  5. the user can force the deletion of endpoints; this is only done
     if condition 3 is met.
  6. when objects are written to a stream we distinguish between
     points and segments; points start by a 'p' and segments by an 's'.
*/


namespace CGAL {

template<class Gt,
	 class ST = Segment_Delaunay_graph_storage_traits_2<Gt>,
	 class D_S = Triangulation_data_structure_2 <
                Segment_Delaunay_graph_vertex_base_2<ST>,
                Segment_Delaunay_graph_face_base_2<Gt> >,
	 class LTag = Tag_false >
class Segment_Delaunay_graph_Linf_2
  : public Segment_Delaunay_graph_2< Gt, ST, D_S, LTag >
{
protected:
  typedef Segment_Delaunay_graph_Linf_2<Gt,ST,D_S,LTag>  Self;
  typedef Segment_Delaunay_graph_2<Gt,ST,D_S,LTag>       Base;

public:
  typedef Gt                                     Geom_traits;
  typedef ST                                     Storage_traits;

  typedef typename Base::Site_2                  Site_2;
  typedef typename Base::Point_2                 Point_2;

  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Base::All_faces_iterator All_faces_iterator;
  typedef typename Base::All_vertices_iterator All_vertices_iterator;
  typedef typename Base::All_edges_iterator All_edges_iterator;

  typedef typename Base::Vertex_handle  Vertex_handle;
  typedef typename Base::Face_handle    Face_handle;
  typedef typename Base::Edge           Edge;
  typedef typename Base::Face           Face;

  typedef typename Base::Vertex_circulator     Vertex_circulator;
  typedef typename Base::Face_circulator       Face_circulator;

protected:
  typedef typename Base::Intersections_tag Intersections_tag;
  typedef typename Base::Face_pair        Face_pair;
  typedef typename Base::Storage_site_2   Storage_site_2;

  typedef typename Base::Vertex_triple    Vertex_triple;


protected:
  using Base::split_storage_site;
  using Base::merge_info;
  using Base::oriented_side;

public:
  using Base::cw;
  using Base::ccw;
  using Base::dimension;
  using Base::number_of_vertices;
  using Base::number_of_faces;
  using Base::number_of_input_sites;
  using Base::number_of_output_sites;
  using Base::all_vertices_begin;
  using Base::all_vertices_end;
  using Base::all_edges_begin;
  using Base::all_edges_end;
  using Base::all_faces_begin;
  using Base::all_faces_end;
  using Base::tds;
  using Base::incident_vertices;
  using Base::incident_faces;
  using Base::infinite_vertex;
  using Base::is_infinite;
  using Base::storage_traits;

private:
  // CREATION helper
  template<class ITag>
  inline
  void setup_if_intersecting_pointer(ITag tag) {
    setup_if_intersecting_pointer_with_tag(tag);
  }

  void setup_if_intersecting_pointer_with_tag(Tag_false) {
    Base::insert_point_on_segment_ptr = NULL;
  }

  void setup_if_intersecting_pointer_with_tag(Tag_true) {
    Base::insert_point_on_segment_ptr =
      static_cast<typename Base::Insert_on_Type>(
        &Self::insert_point_on_segment);
  }

  void setup_insert_on_pointers_linf(void) {
    Intersections_tag itag;
    setup_if_intersecting_pointer(itag);
    Base::insert_exact_point_on_segment_ptr =
      static_cast<typename Base::Insert_Exact_on_Type>(
        &Self::insert_exact_point_on_segment);
  }

public:
  // CREATION
  //---------
  Segment_Delaunay_graph_Linf_2(const Geom_traits& gt = Geom_traits(),
			   const Storage_traits& st = Storage_traits())
    : Base(gt,st)
  {
    setup_insert_on_pointers_linf();
  }

  template< class Input_iterator >
  Segment_Delaunay_graph_Linf_2(Input_iterator first, Input_iterator beyond,
			   const Geom_traits& gt = Geom_traits(),
			   const Storage_traits& st = Storage_traits())
    : Base(gt, st)
  {
    setup_insert_on_pointers_linf();
    Base::insert(first, beyond);
  }

  Segment_Delaunay_graph_Linf_2(const Self& other)
    : Base(other) {}

  Self& operator=(const Self& other)
  {
    if ( this != &other ) {
      Base::operator=(other);
    }
    return (*this);
  }


public:

// print face in standard output
  void face_output(const char *before, Face_handle f,
                   const char *after) const;

  Oriented_side
  oriented_side_face_tiebreak(
    Face_handle f, const Vertex_handle& v, const Site_2 & sitev,
    const Site_2 & sitev_supp, const Site_2 & t) const;


  Face_pair
  find_faces_to_split(const Vertex_handle& v, const Site_2& t,
                      unsigned int & flips_nop,
                      unsigned int & flips_pon) const;

  Vertex_triple
  insert_point_on_segment(const Storage_site_2& ss, const Site_2& t,
			  Vertex_handle v, const Tag_true&);

  Vertex_triple
  insert_exact_point_on_segment(const Storage_site_2& ss, const Site_2& t,
				Vertex_handle v);


  // Linf tiebreaker for finite vertex
  inline Oriented_side
  oriented_side(const Site_2& s1, const Site_2& s2, const Site_2& s3,
		const Site_2& supp, const Site_2& p,
                const Point_2& pt ) const {
    CGAL_precondition( supp.is_segment() && p.is_point() );
    return Base::geom_traits().oriented_side_2_object()(
        s1, s2, s3, supp, p, pt);
  }

  // Linf tiebreaker for infinite vertex
  inline Oriented_side
  oriented_side(const Site_2& s1, const Site_2& s2,
		const Site_2& supp, const Site_2& p,
                const Point_2& pt ) const {
    CGAL_precondition( supp.is_segment() && p.is_point() );
    return Base::geom_traits().oriented_side_2_object()(
        s1, s2, supp, p, pt);
  }


  void file_output_verbose(std::ostream& os) const {
    const char inf_vertex[] = "infinite vertex";
    const char v_id[] = {'p', 'q', 'r', 's'};

    os << "SDG verbose output" << std::endl;
    os << "==================" << std::endl;
    os << "dimension of sdg = " << dimension() << std::endl;
    os << "number of vertices in sdg = "
      << number_of_vertices() << std::endl;
    os << "number of faces in sdg = " << number_of_faces() << std::endl;
    os << "number of input sites in sdg = "
      << number_of_input_sites() << std::endl;
    os << "number of output sites in sdg = "
      << number_of_output_sites() << std::endl;
    // print the vertices of the segment Delaunay graph
    os << std::endl;
    os << "Output sites / vertices of sdg:" << std::endl;
    os << "-------------------------------" << std::endl;
    Finite_vertices_iterator     fvit;
    All_vertices_iterator avit;
    Vertex_handle vh;

    if (dimension() < 2) {
      CGAL_assertion(number_of_faces() == 0);
    }

    int v_count=0;
    for (avit = all_vertices_begin();
         avit != all_vertices_end();
         ++avit) {
      vh = avit;
      if (is_infinite(vh)) {
        os << "vertex " << ++v_count << " : " << inf_vertex << std::endl;
      } else {
        os << "vertex " << ++v_count << " : " << avit->site() << std::endl;
      }
    }

    os << std::endl;
    os << "Edges of sdg:" << std::endl;
    os << "-------------" << std::endl;

    All_edges_iterator eit = all_edges_begin();
    for (int k = 1; eit != all_edges_end(); ++eit, ++k) {
      Edge e = *eit;
      // get the vertices defining the edge

      /*Vertex_handle v[] = {
        (e.first)->vertex(ccw(e.second)),
        (this->dimension() == 1 ?
          (e.first)->vertex(cw(e.second)):
          (e.first)->vertex(e.second))
      };*/
      Vertex_handle v[] = {
        (e.first)->vertex(ccw(e.second)),
        (e.first)->vertex(cw(e.second))
      };


      os << "--- Edge " << k << " ---" << std::endl;
      for (int i = 0; i < 2; i++) {
        // check if the vertex is the vertex at infinity; if yes, print
        // the corresponding string, otherwise print the site
        if ( is_infinite(v[i]) ) {
          os << v_id[i] << ": " << inf_vertex << std::endl;
        } else {
          os << v_id[i] << ": " << v[i]->site() << std::endl;
        }
      }
      os << std::endl;
    }

    os << std::endl;
    os << "Faces of sdg:" << std::endl;
    os << "-------------" << std::endl;

    All_faces_iterator fit = all_faces_begin();
    for (int k = 1; fit != all_faces_end(); ++fit, ++k) {
      Face f = *fit;
      // get the vertices defining the face
      Vertex_handle v[] = { f.vertex(0),
        f.vertex(1),
        f.vertex(2)
      };

      os << "--- Face " << k << " ---" << std::endl;
      for (int i = 0; i < 3; i++) {
        // check if the vertex is the vertex at infinity; if yes, print
        // the corresponding string, otherwise print the site
        if ( is_infinite(v[i]) ) {
          os << v_id[i] << ": " << inf_vertex << std::endl;
        } else {
          os << v_id[i] << ": " << v[i]->site() << std::endl;
        }
      }
      os << std::endl;
    }

    if (dimension() > 1) {
      os << std::endl;
      os << "Voronoi edges of sdg:" << std::endl;
      os << "---------------------" << std::endl;

      All_edges_iterator efit = all_edges_begin();
      for (int k = 1; efit != all_edges_end(); ++efit, ++k) {
        Edge e = *efit;
        // get the vertices defining the Voronoi edge
        Vertex_handle v[] = {
          (e.first)->vertex( ccw(e.second) ),
          (e.first)->vertex(  cw(e.second) ),
          (e.first)->vertex(     e.second  ),
          (this->dimension() == 1) ?
            (e.first)->vertex(   e.second  ) :
            tds().mirror_vertex(e.first, e.second)
        };

        os << "--- Voronoi Edge " << k << " ---" << std::endl;
        for (int i = 0; i < 4; i++) {
          // check if the vertex is the vertex at infinity; if yes, print
          // the corresponding string, otherwise print the site
          if ( is_infinite(v[i]) ) {
            os << v_id[i] << ": " << inf_vertex << std::endl;
          } else {
            os << v_id[i] << ": " << v[i]->site() << std::endl;
          }
        }
        os << std::endl;
      }
    }

    // a counterclockwise traversal of the vertices adjacent to the
    // infinite_vertex is a clockwise traversal of the convex hull
    os << std::endl;
    os << "Convex-hull of sdg:" << std::endl;
    os << "-------------------" << std::endl;
    Vertex_handle v_inf = infinite_vertex();
    Vertex_circulator	 vc1 = incident_vertices(v_inf);
    Vertex_circulator	 vc2 = vc1;
    CGAL_precondition(is_infinite(v_inf));
    /*if (is_infinite(v_inf)){
      os << "vertex 0 : " << inf_vertex << std::endl;
    }*/
    if (vc1 != 0) {
      unsigned int cnt = 0;
      do {
        vh = vc1;
        CGAL_assertion(! is_infinite(vh));
        if (is_infinite(vh)) {
          os << "vertex " << ++cnt << " : " << inf_vertex << std::endl;
        } else {
          os << "vertex " << ++cnt << " : " << vc1->site() << std::endl;
        }
      } while (++vc1 != vc2);
    }

    os << "=======================" << std::endl;
    os << "SDG verbose output ends" << std::endl;
    os << "=======================" << std::endl;

  }

};

} //namespace CGAL


#include <CGAL/Segment_Delaunay_graph_Linf_2/Segment_Delaunay_graph_Linf_2_impl.h>


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_LINF_2_H
