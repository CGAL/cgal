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

//**************************************************************************
//**************************************************************************

template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2(const Geom_traits& traits)
  : Base(traits), random((long)0)
{ 
  hierarchy[0] = this; 
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(traits);
}


// copy constructor duplicates vertices and faces
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
Segment_Voronoi_diagram_hierarchy_2
(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
    : Base(), random((long)0)
{ 
  // create an empty triangulation to be able to delete it !
  hierarchy[0] = this; 
  for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i)
    hierarchy[i] = new Base(svd.geom_traits());
  copy_triangulation(svd);
} 
 

//Assignement
template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
operator=(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
{
  copy_triangulation(svd);
  return *this;
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::   
copy_triangulation
(const Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag> &svd)
{
  std::map< Vertex_handle, Vertex_handle > V;
  {
    for(int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
      *(hierarchy[i]) = *svd.hierarchy[i];
    }
  }

  //up and down have been copied in straightforward way
  // compute a map at lower level
  {
    for(All_vertices_iterator it = hierarchy[0]->all_vertices_begin(); 
	it != hierarchy[0]->all_vertices_end(); ++it) {
      if ( it->up() != Vertex_handle() ) {
	V[ it->up()->down() ] = it;
      }
    }
  }
  {
    for(int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
      for(All_vertices_iterator it = hierarchy[i]->all_vertices_begin(); 
	  it != hierarchy[i]->all_vertices_end(); ++it) {
	// down pointer goes in original instead in copied triangulation
	it->set_down(V[it->down()]);
	// make reverse link
	it->down()->set_up( it );
	// make map for next level
	if ( it->up() != Vertex_handle() ) {
	  V[ it->up()->down() ] = it;
	}
      }
    }
  }

  // copy the point container
  hierarchy[0]->pc_ = svd.hierarchy[0]->pc_;
}

template<class Gt, class STag, class PC, class DS, class LTag>
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
~Segment_Voronoi_diagram_hierarchy_2()
{
  clear();
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i){ 
    delete hierarchy[i];
  }
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
clear()
{
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    hierarchy[i]->clear();
  }
}

//------------------
// INSERTION METHODS
//------------------

template<class Gt, class STag, class PC, class DS, class LTag>
inline typename 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_point(const Point_2& p, int level)

{
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Vertex_handle vertices[svd_hierarchy_2__maxlevel];
  
  insert_point(p, level, vertices);

  return vertices[0];
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_point(const Point_2& p, int level, Vertex_handle* vertices)
{
  CGAL_precondition( level != UNDEFINED_LEVEL );

  Vertex_handle vertex;
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];

  nearest_neighbor(p, vnear, false);

  vertex = hierarchy[0]->insert(p, vnear[0]);

  if ( vertices != NULL ) { vertices[0] = vertex; }

  CGAL_assertion( vertex != Vertex_handle() );

  // insert at other levels
  Vertex_handle previous = vertex;
      
  int k = 1;
  while ( k <= level ) {
    vertex = hierarchy[k]->insert(p, vnear[k]);

    CGAL_assertion( vertex != Vertex_handle() );

    if ( vertices != NULL ) { vertices[k] = vertex; }

    vertex->set_down(previous); // link with other levels
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }
}


template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment_with_tag(const Point_2& p0, const Point_2& p1,
			int level, Tag_false stag)
{
  // the tag is false so we do NOT insert segments in hierarchy
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site_2 t(p0, p1);

  if ( is_degenerate_segment(t) ) {
    return insert_point(p0, level);
  }

  Vertex_handle vertices0[svd_hierarchy_2__maxlevel];
  Vertex_handle vertices1[svd_hierarchy_2__maxlevel];

  insert_point(p0, level, vertices0);
  insert_point(p1, level, vertices1);

  CGAL_assertion( vertices0[0] != Vertex_handle() );
  CGAL_assertion( vertices1[0] != Vertex_handle() );

  Storage_site_2 ss = create_storage_site(vertices0[0], vertices1[0]);

  Vertex_handle vertex;

  if ( hierarchy[0]->number_of_vertices() == 2 ) {
    vertex = hierarchy[0]->insert_third(vertices0[0], vertices1[0]);
  } else {
    vertex = insert_segment_interior(t, ss, vertices0[0], level, stag);
  }

  return vertex;
}



template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment_with_tag(const Point_2& p0, const Point_2& p1,
			int level, Tag_true stag)
{
  // the tag is true so we DO insert segments in hierarchy
  if ( level == UNDEFINED_LEVEL ) {
    level = random_level();
  }

  Site_2 t(p0, p1);

  if ( is_degenerate_segment(t) ) {
    return insert_point(p0, level);
  }

  Vertex_handle vertices0[svd_hierarchy_2__maxlevel];
  Vertex_handle vertices1[svd_hierarchy_2__maxlevel];

  insert_point(p0, level, vertices0);
  insert_point(p1, level, vertices1);

  CGAL_assertion( vertices0[0] != Vertex_handle() );
  CGAL_assertion( vertices1[0] != Vertex_handle() );

  Storage_site_2 ss = create_storage_site(vertices0[0], vertices1[0]);

  Vertex_handle vertex;

  if ( hierarchy[0]->number_of_vertices() == 2 ) {
    vertex = hierarchy[0]->insert_third(vertices0[0], vertices1[0]);

    Vertex_handle previous = vertex;
    Vertex_handle first = vertex;
      
    int k = 1;
    while (k <= level ){
      CGAL_assertion( hierarchy[k]->number_of_vertices() == 2 );

      vertex = hierarchy[k]->insert_third(vertices0[k], vertices1[k]);

      CGAL_assertion( vertex != Vertex_handle() );

      vertex->set_down(previous); // link with level above
      previous->set_up(vertex);
      previous = vertex;
      k++;
    }

    vertex = first;
  } else {
    vertex = insert_segment_interior(t, ss, vertices0, vertices1,
				     level, stag);
  }

  return vertex;
}

//========================================================================

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
			Vertex_handle* vertices0,
			Vertex_handle* vertices1,
			int level, Tag_true stag)
{
  // insert the interior of a segment, and DO insert segments in
  // upper levels of the hierarchy
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() >= 2 );

  CGAL_assertion( vertices0[0] != Vertex_handle() );
  // MK: add here code that checks if the inserted segment has already
  // been inserted; MAYBE THIS IS NOT NEEDED; I ALREADY DO IT IN
  // do_intersect

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vertices0[0]->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( same_segments(t, vv) ) {
      return vv;
    }
    if ( do_intersect(t, vv) ) {
      if ( t.is_segment() ) {
	Intersections_tag itag;
	return insert_intersecting_segment_with_tag(ss, t, vv, level,
						    itag, stag);
      }
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = vertices0[0]->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  std::map<Face_handle,Sign> sign_map;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    sign_map[f] = s;

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // segments must have a conflict with at least one vertex
  CGAL_assertion( s == NEGATIVE );

  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  hierarchy[0]->initialize_conflict_region(start_f, l);
  hierarchy[0]->expand_conflict_region(start_f, t, ss, l, fm,
				       sign_map, vcross, NULL);

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      Intersections_tag itag;
      return insert_intersecting_segment_with_tag(ss, t, vcross.second,
						  level, itag, stag);
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = hierarchy[0]->create_vertex(ss);

  hierarchy[0]->retriangulate_conflict_region(v, l, fm);

  // insert at other levels
  Vertex_handle previous = v;
  Vertex_handle vertex = v;

  int k = 1;
  while (k <= level ){
    if ( hierarchy[k]->number_of_vertices() == 2 ) {
      CGAL_precondition(vertices0 != NULL );
      // MK::ERROR: the if-statement below is a hack. I should deal
      //  with this problem homehow else
      if ( vertices1 == NULL ) {
	Vertex_handle v0(hierarchy[k]->finite_vertices_begin());
	Vertex_handle v1(++(hierarchy[k]->finite_vertices_begin()));
	CGAL_precondition( v0 != Vertex_handle() &&
			   v1 != Vertex_handle() );
	vertex = hierarchy[k]->insert_third(v0, v1);
      } else {
	vertex = hierarchy[k]->insert_third(vertices0[k],
					    vertices1[k]);
      }
    } else {
      vertex = hierarchy[k]->insert_segment_interior(t, ss,
						     vertices0[k], false);
    }

    CGAL_assertion( vertex != Vertex_handle() );

    vertex->set_down(previous); // link with level above
    previous->set_up(vertex);
    previous = vertex;
    k++;
  }

  return v;
}


template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_segment_interior(const Site_2& t, const Storage_site_2& ss,
			Vertex_handle vnearest, int level,
			Tag_false stag)
{
  // insert the interior of a segment, but do not insert segment in
  // upper levels of the hierarchy
  CGAL_precondition( t.is_segment() );
  CGAL_precondition( number_of_vertices() >= 2 );

  CGAL_assertion( vnearest != Vertex_handle() );
  // MK: add here code that checks if the inserted segment has already
  // been inserted; MAYBE THIS IS NOT NEEDED; I ALREADY DO IT IN
  // do_intersect

  // find the first conflict

  // first look if there are intersections...
  Vertex_circulator vc = vnearest->incident_vertices();
  Vertex_circulator vc_start = vc;
  do {
    Vertex_handle vv(vc);
    if ( same_segments(t, vv) ) {
      return vv;
    }
    if ( do_intersect(t, vv) ) {
      if ( t.is_segment() ) {
	Intersections_tag itag;
	return insert_intersecting_segment_with_tag(ss, t, vv, level,
						    itag, stag);
      }
    }
    ++vc;
  } while ( vc != vc_start );

  // first look for conflict with vertex
  Face_circulator fc_start = vnearest->incident_faces();
  Face_circulator fc = fc_start;
  Face_handle start_f;
  Sign s;

  std::map<Face_handle,Sign> sign_map;

  do {
    Face_handle f(fc);

    s = incircle(f, t);

    sign_map[f] = s;

    if ( s == NEGATIVE ) {
      start_f = f;
      break;
    }
    ++fc;
  } while ( fc != fc_start );

  // segments must have a conflict with at least one vertex
  CGAL_assertion( s == NEGATIVE );

  // we are in conflict with a Voronoi vertex; start from that and 
  // find the entire conflict region and then repair the diagram
  List l;
  Face_map fm;

  std::pair<bool, Vertex_handle> vcross(false, Vertex_handle());

  // MK:: NEED TO WRITE A FUNCTION CALLED find_conflict_region WHICH
  // IS GIVEN A STARTING FACE, A LIST, A FACE MAP, A VERTEX MAP AND A
  // LIST OF FLIPPED EDGES AND WHAT IS DOES IS INITIALIZE THE CONFLICT 
  // REGION AND EXPANDS THE CONFLICT REGION.
  hierarchy[0]->initialize_conflict_region(start_f, l);
  hierarchy[0]->expand_conflict_region(start_f, t, ss, l, fm,
				       sign_map, vcross, NULL);

  // the following condition becomes true only if intersecting
  // segments are found
  if ( vcross.first ) {
    if ( t.is_segment() ) {
      Intersections_tag itag;
      return insert_intersecting_segment_with_tag(ss, t, vcross.second,
						  level, itag, stag);
    }
  }

  // no intersecting segment has been found; we insert the segment as
  // usual...
  Vertex_handle v = hierarchy[0]->create_vertex(ss);

  hierarchy[0]->retriangulate_conflict_region(v, l, fm);

  return v;
}

//========================================================================

//--------------------------------------------------------------------
// insertion of intersecting site
//--------------------------------------------------------------------

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int, Tag_false itag, Tag_false stag)
{
  static int i = 0;
  if ( i == 0 ) {
    i = 1;
    std::cerr << std::endl;
    std::cerr << "WARNING:" << std::endl;
    std::cerr << "A segment-segment intersection was found."
	      << std::endl;
    std::cerr << "The segment Voronoi diagram class is not configured"
	      << " to handle this situation." << std::endl;
    std::cerr << "Please look at the documentation on how to handle"
	      << " this behavior." << std::endl;
    std::cerr << std::endl;
  }
  return Vertex_handle();
}


template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int, Tag_false itag, Tag_true stag)
{
  static int i = 0;
  if ( i == 0 ) {
    i = 1;
    std::cerr << std::endl;
    std::cerr << "WARNING:" << std::endl;
    std::cerr << "A segment-segment intersection was found."
	      << std::endl;
    std::cerr << "The segment Voronoi diagram class is not configured"
	      << " to handle this situation." << std::endl;
    std::cerr << "Please look at the documentation on how to handle"
	      << " this behavior." << std::endl;
    std::cerr << std::endl;
  }
  return Vertex_handle();
}


template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_false stag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  if ( same_segments(t, v->site()) ) {
    return v;
  }

  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx( ss.point_handle(0), ss.point_handle(1),
		      ssitev.point_handle(0), ssitev.point_handle(1) );

  Site_2 sitev(v->site());
  Site_2 sx(t.point(0), t.point(1), sitev.point(0), sitev.point(1));

  Face_circulator fc1 = incident_faces(v);
  Face_circulator fc2 = fc1; ++fc2;
  Face_circulator fc_start = fc1;
  Face_handle f1, f2;
  bool found_f1 = false, found_f2 = false;
  Site_2 sitev_supp(sitev.point(0), sitev.point(1));
  do {
    Face_handle ff1(fc1), ff2(fc2);
    Oriented_side os1 = oriented_side(fc1->vertex(0)->site(),
				      fc1->vertex(1)->site(),
				      fc1->vertex(2)->site(),
				      sitev_supp, sx);
    Oriented_side os2 = oriented_side(fc2->vertex(0)->site(),
				      fc2->vertex(1)->site(),
				      fc2->vertex(2)->site(),
				      sitev_supp, sx);
    if ( !found_f1 &&
	 os1 != ON_POSITIVE_SIDE && os2 == ON_POSITIVE_SIDE ) {
      f1 = ff2;
      found_f1 = true;
    }

    if ( !found_f2 &&
	 os1 == ON_POSITIVE_SIDE && os2 != ON_POSITIVE_SIDE ) {
      f2 = ff2;
      found_f2 = true;
    }

    if ( found_f1 && found_f2 ) { break; }

    ++fc1, ++fc2;
  } while ( fc_start != fc1 ); 

  CGAL_assertion( f1 != f2 );

  Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
    qq = hierarchy[0]->_tds.split_vertex(v, f1, f2);

  // now I need to update the sites for vertices v1 and v2
  Vertex_handle v1 = qq.first;
  Storage_site_2 ssv1;
  Site_2 sv1;
  if ( sitev.is_exact(0) ) {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), true);
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), true);
  } else {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    sitev.point(2), sitev.point(3),
		    t.point(0), t.point(1));
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ssitev.point_handle(2), ssitev.point_handle(3),
		     ss.point_handle(0), ss.point_handle(1));
  }
  v1->set_site( ssv1 );

  Vertex_handle v2 = qq.second;
  Storage_site_2 ssv2;
  Site_2 sv2;
  if ( sitev.is_exact(1) ) {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), false);
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), false);
  } else {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1),
		    sitev.point(4), sitev.point(5));
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1),
		     ssitev.point_handle(4), ssitev.point_handle(5));
  }
  v2->set_site( ssv2 );

  Vertex_handle vsx =
    hierarchy[0]->_tds.insert_in_edge(qq.third, cw(qq.third->index(v1)));

  vsx->set_site(ssx);

  //  insert_point(sx, level, vertices);

  // use insert_point(sx, level, vertices) instead, but this version
  // of insert_point must be enriched by adding not only the max
  // level, but also the min level...
  // MK::ERROR: I may need to do geometric filtering in the
  // Orientation_2 and Incircle_2 predicates for points, because now I
  // insert lots of points that are collinear...
  {
    Vertex_handle vnear[svd_hierarchy_2__maxlevel];

    nearest_neighbor(sx, vnear, false);

    Vertex_handle vertex = vsx;

    Vertex_handle previous = vertex;

    // MK::ERROR: Another idea (since I am not inserting any segments
    // in the first place, is to recompute the level of each of these
    // points of intersection. This way I do not get long "lines" of
    // points appearing in the upper levels...; this what is done with
    // the new_level variable below...
    int new_level = random_level();
    int k = 1;
    while ( k <= new_level ) {
      // MK::ERROR: this is a hack. I need to change the code in the
      // segment Voronoi diagram class, so that I can insert sites as
      // the first, second, or third site...; actually this is
      // problematic only if the number of vertices is exactly 2; it
      // cannot be smaller since we already have added the endpoints
      // of the segment at level k.
      if ( hierarchy[k]->number_of_vertices() <= 2 ) {
	//	CGAL_assertion( hierarchy[k]->number_of_vertices() == 2 );
	break;
      }
      vertex = hierarchy[k]->insert_point(ssx, sx, vnear[k]);

      CGAL_assertion( vertex != Vertex_handle() );

      vertex->set_down(previous); // link with other levels
      previous->set_up(vertex);
      previous = vertex;
      k++;
    }
  }

  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_exact(0) ) {
    s3.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), true);
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1), true);
  } else {
    s3.set_segment(t.point(0), t.point(1),
		   t.point(2), t.point(3),
		   sitev.point(0), sitev.point(1));
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ss.point_handle(2), ss.point_handle(3),
		    ssitev.point_handle(0), ssitev.point_handle(1));
  }

  if ( t.is_exact(1) ) {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), false);
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1), false);
  } else {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1),
		   t.point(4), t.point(5));
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1),
		    ss.point_handle(4), ss.point_handle(5));
  }

  insert_segment_interior(s3, ss3, vsx, level, stag);
  insert_segment_interior(s4, ss4, vsx, level, stag);
  return vsx;
}

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				     const Site_2& t, Vertex_handle v,
				     int level,
				     Tag_true itag, Tag_true stag)
{
  CGAL_precondition( t.is_segment() && v->is_segment() );

  // MK::ERROR: I have to remove this; too expensive...
  CGAL_precondition( do_intersect(t, v->site()) );

  if ( same_segments(t, v->site()) ) {
    // MK::ERROR: I may need to insert it to levels higher than its
    // previous level...
    return v;
  }

  Storage_site_2 ssitev = v->storage_site();
  Storage_site_2 ssx( ss.point_handle(0), ss.point_handle(1),
		      ssitev.point_handle(0), ssitev.point_handle(1) );

  Site_2 sitev(v->site());
  Site_2 sx(t.point(0), t.point(1), sitev.point(0), sitev.point(1));
  Site_2 sitev_supp(sitev.point(0), sitev.point(1));

  Storage_site_2 ssv1;
  Site_2 sv1;
  if ( sitev.is_exact(0) ) {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), true);
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), true);
  } else {
    sv1.set_segment(sitev.point(0), sitev.point(1),
		    sitev.point(2), sitev.point(3),
		    t.point(0), t.point(1));
    ssv1.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ssitev.point_handle(2), ssitev.point_handle(3),
		     ss.point_handle(0), ss.point_handle(1));
  }

  Storage_site_2 ssv2;
  Site_2 sv2;
  if ( sitev.is_exact(1) ) {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1), false);
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1), false);
  } else {
    sv2.set_segment(sitev.point(0), sitev.point(1),
		    t.point(0), t.point(1),
		    sitev.point(4), sitev.point(5));
    ssv2.set_segment(ssitev.point_handle(0), ssitev.point_handle(1),
		     ss.point_handle(0), ss.point_handle(1),
		     ssitev.point_handle(4), ssitev.point_handle(5));
  }

  Vertex_handle vertices[svd_hierarchy_2__maxlevel];

  int v_level = find_level(v);

  Vertex_handle vertex = v;
  Vertex_handle v1_old, v2_old, vsx_old;

  int m = 0;
  while ( m <= v_level ) {
    // MK::ERROR: I have to remove this; too expensive...
    CGAL_precondition( do_intersect(t, vertex->site()) );

    Face_circulator fc1 = hierarchy[m]->incident_faces(vertex);
    Face_circulator fc2 = fc1; ++fc2;
    Face_circulator fc_start = fc1;
    Face_handle f1, f2;
    bool found_f1 = false, found_f2 = false;
    do {
      Face_handle ff1(fc1), ff2(fc2);
      CGAL_assertion( !is_infinite(ff1) && !is_infinite(ff2) );
      Oriented_side os1 = oriented_side(fc1->vertex(0)->site(),
					fc1->vertex(1)->site(),
					fc1->vertex(2)->site(),
					sitev_supp, sx);
      Oriented_side os2 = oriented_side(fc2->vertex(0)->site(),
					fc2->vertex(1)->site(),
					fc2->vertex(2)->site(),
					sitev_supp, sx);
      if ( !found_f1 &&
	   os1 != ON_POSITIVE_SIDE && os2 == ON_POSITIVE_SIDE ) {
	f1 = ff2;
	found_f1 = true;
      }

      if ( !found_f2 &&
	   os1 == ON_POSITIVE_SIDE && os2 != ON_POSITIVE_SIDE ) {
	f2 = ff2;
	found_f2 = true;
      }

      if ( found_f1 && found_f2 ) { break; }

      ++fc1, ++fc2;
    } while ( fc_start != fc1 ); 
    
    CGAL_assertion( f1 != f2 );
    CGAL_assertion( !is_infinite(f1) && !is_infinite(f2) );

    Vertex_handle vertex_up = vertex->up();

    Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
      qq = hierarchy[m]->_tds.split_vertex(vertex, f1, f2);

    // now I need to update the sites for vertices v1 and v2
    Vertex_handle v1 = qq.first;
    Vertex_handle v2 = qq.second;
    v1->set_site( ssv1 );
    v2->set_site( ssv2 );

    CGAL_assertion( v1->is_segment() && v2->is_segment() );

    Vertex_handle vsx =
      hierarchy[m]->_tds.insert_in_edge(qq.third,
					cw(qq.third->index(v1)));
    vsx->set_site(ssx);

    if ( m > 0 ) {
      if ( same_segments(v1->site(), v1_old->site()) ) {
	v1->set_down(v1_old);
	v2->set_down(v2_old);
	v1_old->set_up(v1);
	v2_old->set_up(v2);
      } else {
	v1->set_down(v2_old);
	v2->set_down(v1_old);
	v1_old->set_up(v2);
	v2_old->set_up(v1);
      }
      vsx_old->set_up(vsx);
      vsx->set_down(vsx_old);
    }

    v1_old = v1;
    v2_old = v2;
    vsx_old = vsx;

    vertices[m] = vsx;

    vertex = vertex_up;
    m++;
  }

  //  insert_point(sx, level, vertices);

  // use insert_point(sx, level, vertices) instead, but this version
  // of insert_point must be enriched by adding not only the max
  // level, but also the min level...
  // MK::ERROR: I may need to do geometric filtering in the
  // Orientation_2 and Incircle_2 predicates for points, because now I
  // insert lots of points that are collinear...
  if ( v_level < level ) {

    Vertex_handle vnear[svd_hierarchy_2__maxlevel];

    nearest_neighbor(sx, vnear, false);

    Vertex_handle previous = vertices[v_level];

    // MK::ERROR: Another idea (since I am not inserting any segments
    // in the first place, is to recompute the level of each of these
    // points of intersection. This way I do not get long "lines" of
    // points appearing in the upper levels...; this what is done with
    // the new_level variable below...
    int k = v_level + 1;
    while ( k <= level ) {
      // MK::ERROR: this is a hack. I need to change the code in the
      // segment Voronoi diagram class, so that I can insert sites as
      // the first, second, or third site...; actually this is
      // problematic only if the number of vertices is exactly 2; it
      // cannot be smaller since we already have added the endpoints
      // of the segment at level k.
      if ( hierarchy[k]->number_of_vertices() <= 2 ) {
	CGAL_assertion( hierarchy[k]->number_of_vertices() == 2 );
	// MK::ERROR: this creates a site that is not exact (I am
	// approximating an intersection point site by a regular point
	// site...
	//	CGAL_assertion( false );
	vertex = hierarchy[k]->insert_third(sx, ssx);
      } else {
	vertex = hierarchy[k]->insert_point(ssx, sx, vnear[k]);
      }

      vertices[k] = vertex;

      CGAL_assertion( vertex != Vertex_handle() );

      vertex->set_down(previous); // link with other levels
      previous->set_up(vertex);
      previous = vertex;
      k++;
    }
  }

  Storage_site_2 ss3, ss4;
  Site_2 s3, s4;
  if ( t.is_exact(0) ) {
    s3.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), true);
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1), true);
  } else {
    s3.set_segment(t.point(0), t.point(1),
		   t.point(2), t.point(3),
		   sitev.point(0), sitev.point(1));
    ss3.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ss.point_handle(2), ss.point_handle(3),
		    ssitev.point_handle(0), ssitev.point_handle(1));
  }

  if ( t.is_exact(1) ) {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1), false);
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1), false);
  } else {
    s4.set_segment(t.point(0), t.point(1),
		   sitev.point(0), sitev.point(1),
		   t.point(4), t.point(5));
    ss4.set_segment(ss.point_handle(0), ss.point_handle(1),
		    ssitev.point_handle(0), ssitev.point_handle(1),
		    ss.point_handle(4), ss.point_handle(5));
  }

  insert_segment_interior(s3, ss3, vertices, NULL, level, stag);
  insert_segment_interior(s4, ss4, vertices, vertices, level, stag);
  return vertices[0];
}


//===========================================================================

template<class Gt, class STag, class PC, class DS, class LTag>
typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::Vertex_handle 
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
nearest_neighbor(const Point_2& p, bool force_point) const
{
  Vertex_handle vnear[svd_hierarchy_2__maxlevel];
  nearest_neighbor(Site_2(p), vnear, force_point);
  return vnear[0];
}

template<class Gt, class STag, class PC, class DS, class LTag>
void
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
nearest_neighbor(const Site_2& p,
		 Vertex_handle vnear[svd_hierarchy_2__maxlevel],
		 bool force_point) const
{
  CGAL_precondition( p.is_point() );

  Vertex_handle nearest;
  int level  = svd_hierarchy_2__maxlevel;

  // find the highest level with enough vertices
  while ( hierarchy[--level]->number_of_vertices() 
	  < svd_hierarchy_2__minsize ) {
    if ( !level ) break;  // do not go below 0
  }
  for (unsigned int i = level + 1; i < svd_hierarchy_2__maxlevel; i++) {
    vnear[i] = Vertex_handle();
  }

  while ( level > 0 ) {
    vnear[level] = nearest =
      hierarchy[level]->nearest_neighbor(p, nearest);  

    CGAL_assertion( !hierarchy[level]->is_infinite(vnear[level]) );
    CGAL_assertion( vnear[level] != Vertex_handle() );
    // go at the same vertex on level below
    nearest = nearest->down();
    --level;
  }
  vnear[0] = hierarchy[level]->nearest_neighbor(p, nearest);
  // at level 0
}

template<class Gt, class STag, class PC, class DS, class LTag>
int
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
random_level()
{
  unsigned int l = 0;
  while ( true ) {
    if ( random(svd_hierarchy_2__ratio) ) break;
    ++l;
  }
  if (l >= svd_hierarchy_2__maxlevel)
    l = svd_hierarchy_2__maxlevel -1;
  return l;
}

template<class Gt, class STag, class PC, class DS, class LTag>
inline typename
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::size_type
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>::
find_level(Vertex_handle v) const
{
  CGAL_precondition( v != Vertex_handle() );
  size_type level = 0;
  Vertex_handle vertex = v;
  while ( vertex->up() != Vertex_handle() ) {
    vertex = vertex->up();
    level++;
  }

  return level;
}

template<class Gt, class STag, class PC, class DS, class LTag>
bool
Segment_Voronoi_diagram_hierarchy_2<Gt,STag,PC,DS,LTag>:: 
is_valid(bool verbose, int level) const
{
  bool result(true);

  //verify correctness of triangulation at all levels
  for(unsigned int i = 0; i < svd_hierarchy_2__maxlevel; ++i) {
    if ( verbose ) {
      std::cout << "Level " << i << ": " << std::flush;
    }
    result = result && hierarchy[i]->is_valid(verbose, level);
    if ( verbose ) {
      std::cout << std::endl;
    }
  }
  //verify that lower level has no down pointers
  for( Finite_vertices_iterator it = hierarchy[0]->finite_vertices_begin(); 
       it != hierarchy[0]->finite_vertices_end(); ++it) {
    result = result && ( it->down() == 0 );
  }

  //verify that other levels has down pointer and reciprocal link is fine
  for(unsigned int i = 1; i < svd_hierarchy_2__maxlevel; ++i) {
    for( Finite_vertices_iterator it = hierarchy[i]->finite_vertices_begin(); 
	 it != hierarchy[i]->finite_vertices_end(); ++it) {
      Vertex_handle vit(it);
      result = result && ( it->down()->up() == vit );
    }
  }
  return result;
}



CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_HIERARCHY_2_H

