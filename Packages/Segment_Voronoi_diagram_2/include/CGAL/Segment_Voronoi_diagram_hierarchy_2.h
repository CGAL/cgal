// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

#include <map>

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <CGAL/Random.h>
#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>


CGAL_BEGIN_NAMESPACE

//--------------------------------------------------------------------
//--------------------------------------------------------------------

// parameterization of the hierarchy
#ifdef CGAL_SVD_HIERARCHY_DEMO
const unsigned int svd_hierarchy_2__ratio    = 3;
const unsigned int svd_hierarchy_2__minsize  = 5;
#else
const unsigned int svd_hierarchy_2__ratio    = 30;
const unsigned int svd_hierarchy_2__minsize  = 20;
#endif
const unsigned int svd_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class Gt, class STag = Tag_false,
	   class PC = std::list<typename Gt::Point_2>,
	   class DS = Segment_Voronoi_diagram_data_structure_2<
              Segment_Voronoi_diagram_hierarchy_vertex_base_2<
                 Segment_Voronoi_diagram_vertex_base_2<Gt,PC,
			     typename Gt::Intersections_tag> >,
              Triangulation_face_base_2<Gt> >,
	   class LTag = Tag_false>
class Segment_Voronoi_diagram_hierarchy_2
  : public Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>
{
public:
  typedef Segment_Voronoi_diagram_2<Gt,PC,DS,LTag>  Base;

  typedef typename Base::Geom_traits        Geom_traits;

  typedef typename Geom_traits::Point_2     Point_2;
  typedef typename Geom_traits::Site_2      Site_2;

  typedef typename Base::Vertex_handle      Vertex_handle;
  typedef typename Base::Face_handle        Face_handle;
  typedef typename Base::Edge               Edge;

  typedef typename Base::Vertex_circulator         Vertex_circulator;
  typedef typename Base::Edge_circulator           Edge_circulator;
  typedef typename Base::Face_circulator           Face_circulator;

  typedef typename Base::All_faces_iterator        All_faces_iterator;
  typedef typename Base::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename Base::All_vertices_iterator     All_vertices_iterator;
  typedef typename Base::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename Base::All_edges_iterator        All_edges_iterator;
  typedef typename Base::Finite_edges_iterator     Finite_edges_iterator;

  typedef typename Base::Point_container           Point_container;
  typedef typename Base::size_type                 size_type;

  typedef typename Base::Intersections_tag         Intersections_tag;

  typedef STag                            Insert_segments_in_hierarchy_tag;

private:
  struct Vertex_iterator {};

  static const int UNDEFINED_LEVEL;

  typedef typename Base::Storage_site_2            Storage_site_2;
  typedef typename Base::List                      List;
  typedef typename Base::Face_map                  Face_map;

private:
  // here is the stack of triangulations which form the hierarchy
  Base*   hierarchy[svd_hierarchy_2__maxlevel];
  Random random; // random generator

public:
  Segment_Voronoi_diagram_hierarchy_2
  (const Geom_traits& traits = Geom_traits());
  Segment_Voronoi_diagram_hierarchy_2
  (const Segment_Voronoi_diagram_hierarchy_2& svd);

  Segment_Voronoi_diagram_hierarchy_2 &operator=
  (const  Segment_Voronoi_diagram_hierarchy_2& svd);

  ~Segment_Voronoi_diagram_hierarchy_2();

  // Helping
  void copy_triangulation
  (const Segment_Voronoi_diagram_hierarchy_2 &svd);

  void clear();

  // CHECKING
  bool is_valid(bool verbose = true, int level = 1) const;

  const Base& diagram(unsigned int i) const  {
    CGAL_precondition( i < svd_hierarchy_2__maxlevel );
    return *hierarchy[i];
  }

public:
   // insertion of a point/segment
  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond) {
    return insert_with_tag(first, beyond, Tag_false());
  }

protected:
  template<class Input_iterator>
  size_type insert_with_tag(Input_iterator first,
			    Input_iterator beyond,
			    Tag_true)
  {
    // MK::ERROR: this changes the data I would have to copy them to
    //            a vector first and then do the suffling thing...
    std::random_shuffle(first, beyond);
    return insert_with_tag(first, beyond, Tag_false());
  }

  template<class Input_iterator>
  size_type insert_with_tag(Input_iterator first,
			    Input_iterator beyond,
			    Tag_false)
  {
    // do it the obvious way: insert them as they come;
    // one might think though that it might be better to first insert
    // all end points and then all segments, or a variation of that.

    size_type n_before = number_of_vertices();
    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
    }
    size_type n_after = number_of_vertices();
    return n_after - n_before;
  }

public:
  template<class Input_iterator, class True_false_tag>
  size_type insert(Input_iterator first, Input_iterator beyond,
		   True_false_tag tag)
  {
    return insert_with_tag(first, beyond, tag);
  }

  Vertex_handle  insert(const Point_2& p) {
    return insert_point(p, UNDEFINED_LEVEL);
  }

  Vertex_handle  insert(const Point_2& p0, const Point_2& p1) {
    return insert_segment_with_tag(p0, p1, UNDEFINED_LEVEL,
				   Insert_segments_in_hierarchy_tag());
  }

  Vertex_handle insert(const Point_2& p, Vertex_handle) {
    return insert(p);
  }

  Vertex_handle insert(const Point_2& p1, const Point_2& p2,
		       Vertex_handle) {
    return insert(p1, p2);
  }

  Vertex_handle  insert(const Site_2& t) {
    CGAL_precondition( t.is_exact() );
    if ( t.is_segment() ) {
      Insert_segments_in_hierarchy_tag stag;
      return insert_segment_with_tag(t.source(), t.target(),
				     UNDEFINED_LEVEL, stag);
    } else if ( t.is_point() ) {
      return insert_point(t.point(), UNDEFINED_LEVEL);
    } else {
      CGAL_precondition ( t.is_defined() );
      return Vertex_handle(); // to avoid compiler error
    }
  }


private:
  Vertex_handle insert_point(const Point_2& p, int level);
  void          insert_point(const Point_2& p, int level,
			     Vertex_handle* vertices);
#if 0
  // not implemented yet
  void          insert_point(const Site_2& p, int level,
			     Vertex_handle* vertices);
#endif

  Vertex_handle insert_segment_with_tag(const Point_2& p0,
					const Point_2& p1,
					int level, Tag_true stag); 

  Vertex_handle insert_segment_with_tag(const Point_2& p0,
					const Point_2& p1,
					int level, Tag_false stag); 


  Vertex_handle insert_segment_interior(const Site_2& t,
					const Storage_site_2& ss,
					Vertex_handle* vertices0,
					Vertex_handle* vertices1,
					int level, Tag_true stag);

  Vertex_handle insert_segment_interior(const Site_2& t,
					const Storage_site_2& ss,
					Vertex_handle vnear,
					int level, Tag_false stag);

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v,
				       int level,
				       Tag_false itag, Tag_false stag);

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v,
				       int level,
				       Tag_false itag, Tag_true stag);

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v,
				       int level,
				       Tag_true itag, Tag_false stag);

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v,
				       int level,
				       Tag_true itag, Tag_true stag);


public:
  // nearest neighbor
  Vertex_handle  nearest_neighbor(const Point_2& p,
				  bool force_point = false) const;

  Vertex_handle  nearest_neighbor(const Point_2& p, Vertex_handle)
  {
    return nearest_neighbor(p);
  }

private:
  void nearest_neighbor(const Site_2& p,
			Vertex_handle vnear[svd_hierarchy_2__maxlevel],
			bool force_point) const; 
  int random_level();

  size_type find_level(Vertex_handle v) const;
};

template<class Gt, class STag, class PC, class DS, class LTag>
const int
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::UNDEFINED_LEVEL = -1;


CGAL_END_NAMESPACE



#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Segment_Voronoi_diagram_hierarchy_2.C>
#endif



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

