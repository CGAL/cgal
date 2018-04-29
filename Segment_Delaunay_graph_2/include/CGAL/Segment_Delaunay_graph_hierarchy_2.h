// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>

#include <CGAL/disable_warnings.h>

#include <map>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Segment_Delaunay_graph_vertex_base_2.h>
#include <CGAL/Segment_Delaunay_graph_hierarchy_vertex_base_2.h>
#include <CGAL/Segment_Delaunay_graph_face_base_2.h>


namespace CGAL {

//--------------------------------------------------------------------
//--------------------------------------------------------------------

// parameterization of the hierarchy
#ifdef CGAL_SDG_HIERARCHY_DEMO
const unsigned int sdg_hierarchy_2__ratio    = 3;
const unsigned int sdg_hierarchy_2__minsize  = 5;
#else
const unsigned int sdg_hierarchy_2__ratio    = 30;
const unsigned int sdg_hierarchy_2__minsize  = 20;
#endif
const unsigned int sdg_hierarchy_2__maxlevel = 5;
// maximal number of points is 30^5 = 24 millions !

//--------------------------------------------------------------------
//--------------------------------------------------------------------

template < class Gt,
	   class ST = Segment_Delaunay_graph_storage_traits_2<Gt>,
	   class STag = Tag_false,
	   class D_S = Triangulation_data_structure_2<
              Segment_Delaunay_graph_hierarchy_vertex_base_2<
		Segment_Delaunay_graph_vertex_base_2<ST> >,
              Segment_Delaunay_graph_face_base_2<Gt> >,
	   class LTag = Tag_false,
           class SDGLx = Segment_Delaunay_graph_2<Gt,ST,D_S,LTag> >
class Segment_Delaunay_graph_hierarchy_2
  : public SDGLx
{
protected:
  typedef Segment_Delaunay_graph_hierarchy_2<Gt,ST,STag,D_S,LTag,SDGLx>  Self;

public:
  // PUBLIC TYPES
  //-------------
  typedef SDGLx    Base;

  typedef typename Base::Geom_traits        Geom_traits;
  typedef typename Base::Storage_traits     Storage_traits;

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

  typedef typename Base::Input_sites_iterator      Input_sites_iterator;
  typedef typename Base::Output_sites_iterator     Output_sites_iterator;

  typedef typename Base::Point_handle              Point_handle;

protected:
  typedef typename Base::Handle_map                Handle_map;
  typedef typename Base::Point_handle_pair         Point_handle_pair;

  using Base::merge_info;
  using Base::same_segments;
  using Base::is_degenerate_segment;
  using Base::convert_info;
  using Base::second_endpoint_of_segment;
  using Base::split_storage_site;
  using Base::first_endpoint_of_segment;
  using Base::incircle;

public:
  using Base::is_infinite;

  typedef typename Base::Point_container           Point_container;
  typedef typename Base::size_type                 size_type;

  typedef typename Base::Intersections_tag         Intersections_tag;

  typedef STag                            Insert_segments_in_hierarchy_tag;
  typedef STag                            Segments_in_hierarchy_tag;

protected:
  // LOCAL TYPES
  //------------
  typedef typename Base::Storage_site_2            Storage_site_2;
  typedef typename Base::List                      List;
  typedef typename Base::Face_map                  Face_map;
  typedef typename Base::Vertex_triple             Vertex_triple;

  typedef typename Base::Arrangement_type          Arrangement_type;
  typedef typename Base::AT2                       AT2;

  enum { UNDEFINED_LEVEL = -1 };

  // here is the stack of triangulations which form the hierarchy
  Base*   hierarchy[sdg_hierarchy_2__maxlevel];
  boost::rand48  random; // random generator
public:
  // CONSTRUCTORS
  //-------------
  Segment_Delaunay_graph_hierarchy_2(const Gt& gt = Gt());

  template<class Input_iterator>
  Segment_Delaunay_graph_hierarchy_2(Input_iterator first,
				     Input_iterator beyond,
				     const Gt& gt=Gt())
    : Base(gt)
  {
    init_hierarchy(gt);
    insert(first, beyond);
  }

  Segment_Delaunay_graph_hierarchy_2(const Self& sdg);
  Self& operator=(const Self& sdg);

  // DESTRUCTOR
  //-----------
  ~Segment_Delaunay_graph_hierarchy_2();

public:
  // ACCESS METHODS
  //---------------
  const Base& diagram(unsigned int i) const  {
    CGAL_precondition( i < sdg_hierarchy_2__maxlevel );
    return *hierarchy[i];
  }

public:
  // INSERTION METHODS
  //------------------
  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond) {
    return insert(first, beyond, Tag_false());
  }

  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond,	Tag_true)
  {
    std::vector<Site_2> site_vec;
    for (Input_iterator it = first; it != beyond; ++it) {
      site_vec.push_back(Site_2(*it));
    }

    typedef std::iterator_traits<Input_iterator> Iterator_traits;
    typedef typename Iterator_traits::difference_type Diff_t;

    boost::random_number_generator<boost::rand48, Diff_t> rng(random);
    CGAL::cpp98::random_shuffle(site_vec.begin(), site_vec.end(),rng);
    return insert(site_vec.begin(), site_vec.end(), Tag_false());
  }

  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond,	Tag_false)
  {
    // do it the obvious way: insert them as they come;
    // one might think though that it might be better to first insert
    // all end points and then all segments, or a variation of that.
    size_type n_before = this->number_of_vertices();
    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
    }
    size_type n_after = this->number_of_vertices();
    return n_after - n_before;
  }

  Vertex_handle  insert(const Point_2& p) {
    Point_handle ph = this->register_input_site(p);
    Storage_site_2 ss = 
      this->st_.construct_storage_site_2_object()(ph);
    return insert_point(p, ss, UNDEFINED_LEVEL);
  }

  Vertex_handle  insert(const Point_2& p0, const Point_2& p1) {
    Point_handle_pair php = this->register_input_site(p0,p1);
    Storage_site_2 ss =
      this->st_.construct_storage_site_2_object()(php.first, php.second);
    Vertex_handle v = insert_segment(p0, p1, ss, UNDEFINED_LEVEL);
    if ( v == Vertex_handle() ) {
      this->unregister_input_site( php.first, php.second );
    }
    return v;
  }

  Vertex_handle insert(const Vertex_handle& v0, const Vertex_handle& v1) {
    return hierarchy[0]->insert(v0, v1);
  }

  Vertex_handle insert(const Point_2& p, Vertex_handle) {
    return insert(p);
  }

  Vertex_handle insert(const Point_2& p0, const Point_2& p1,
		       Vertex_handle) {
    return insert(p0, p1);
  }

  Vertex_handle  insert(const Site_2& t) {
    // the intended use is to unify the calls to insert(...);
    // thus the site must be an exact one; 
    CGAL_precondition( t.is_input() );

    if ( t.is_segment() ) {
      Point_handle_pair php =
	this->register_input_site(t.source(), t.target());
      Storage_site_2 ss =
	this->st_.construct_storage_site_2_object()(php.first, php.second);
      Vertex_handle v =
	insert_segment(t.source(), t.target(), ss, UNDEFINED_LEVEL);
      if ( v == Vertex_handle() ) {
	this->unregister_input_site( php.first, php.second );
      }
      return v;
    } else if ( t.is_point() ) {
      Point_handle ph = this->register_input_site( t.point() );
      Storage_site_2 ss = this->st_.construct_storage_site_2_object()(ph);
      return insert_point(t.point(), ss, UNDEFINED_LEVEL);
    } else {
      CGAL_precondition ( t.is_defined() );
      return Vertex_handle(); // to avoid compiler error
    }
  }

  inline Vertex_handle insert(const Site_2& t, Vertex_handle) {
    return insert(t);
  }

  template<class Info_t>
  inline
  Vertex_handle insert(const Site_2& t, const Info_t& info)
  {
    typedef typename Storage_traits::Info Info;
    CGAL_SEGMENT_DELAUNAY_GRAPH_2_NS::Internal::
      Check_type_equality_for_info<Info_t, Info>();
    // the intended use is to unify the calls to insert(...);
    // thus the site must be an exact one; 
    CGAL_precondition( t.is_input() );

    if ( t.is_segment() ) {
      Point_handle_pair php =
	this->register_input_site(t.source(), t.target());
      Storage_site_2 ss =
	this->st_.construct_storage_site_2_object()(php.first, php.second);
      ss.set_info(info);
      Vertex_handle v =
	insert_segment(t.source(), t.target(), ss, UNDEFINED_LEVEL);
      if ( v == Vertex_handle() ) {
	this->unregister_input_site( php.first, php.second );
      }
      return v;
    } else if ( t.is_point() ) {
      Point_handle ph = this->register_input_site( t.point() );
      Storage_site_2 ss = this->st_.construct_storage_site_2_object()(ph);
      ss.set_info(info);
      return insert_point(t.point(), ss, UNDEFINED_LEVEL);
    } else {
      CGAL_precondition ( t.is_defined() );
      return Vertex_handle(); // to avoid compiler error
    }
  }

  template<class Info_t>
  inline
  Vertex_handle insert(const Site_2& t, const Info_t& info, Vertex_handle)
  {
    return insert(t, info);
  }

protected:
  Vertex_handle insert_point(const Point_2& p, const Storage_site_2& ss,
			     int level) {
    if ( level == UNDEFINED_LEVEL ) {
      level = random_level();
    }

    Vertex_handle vertices[sdg_hierarchy_2__maxlevel];
  
    insert_point(p, ss, level, vertices);

    return vertices[0];
  }

  //  std::pair<bool,Vertex_triple>
  std::pair<bool,int>
                insert_point(const Point_2& p, const Storage_site_2& ss,
			     int level, Vertex_handle* vertices);

  void          insert_point(const Site_2& t, const Storage_site_2& ss,
			     int low, int high, Vertex_handle vbelow,
			     Vertex_handle* vertices);

  Vertex_handle insert_segment(const Point_2& p0, const Point_2& p1,
			       const Storage_site_2& ss, int level); 

  Vertex_handle insert_segment_interior(const Site_2& t,
					const Storage_site_2& ss,
					const Vertex_handle* vertices0,
					int level);

  void insert_segment_in_upper_levels(const Site_2& t,
				      const Storage_site_2& ss,
				      Vertex_handle vbelow,
				      const Vertex_handle* vertices0,
				      int level, Tag_true);

  void insert_segment_in_upper_levels(const Site_2& ,
				      const Storage_site_2& ,
				      Vertex_handle ,
				      const Vertex_handle* ,
				      int , Tag_false) {}

  Vertex_handle insert_segment_on_point(const Storage_site_2& ss,
					const Vertex_handle& v,
					int level, int which);

  template<class Tag>
  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2&,
				       const Site_2& ,
				       Vertex_handle ,
				       int , Tag_false /* itag */, Tag) {
    print_error_message();
    return Vertex_handle();
  }

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
  // REMOVAL METHODS
  //----------------
  bool remove(const Vertex_handle& v);

public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point_2& p,
				  bool force_point = false) const;

  Vertex_handle  nearest_neighbor(const Point_2& p, Vertex_handle) const
  {
    return nearest_neighbor(p);
  }

protected:
  void nearest_neighbor(const Site_2& p,
			Vertex_handle vnear[sdg_hierarchy_2__maxlevel],
			bool force_point) const; 

public:
  // MISCELLANEOUS
  //--------------
  void init_hierarchy(const Geom_traits& gt);

  void copy(const Segment_Delaunay_graph_hierarchy_2& sdgh);

  void swap(Segment_Delaunay_graph_hierarchy_2& sdgh);
  void clear();

public:
  // FILE I/O
  //---------
  void file_input(std::istream&);
  void file_output(std::ostream&) const;

public:
  // VALIDITY CHECK
  //---------------
  bool is_valid(bool verbose = false, int level = 1) const;

protected:
  // LOCAL HELPER METHODS
  //---------------------
  int random_level() {
    boost::geometric_distribution<> proba(1.0/sdg_hierarchy_2__ratio);
    boost::variate_generator<boost::rand48&, boost::geometric_distribution<> > die(random, proba);

    return (std::min)(die(), (int)sdg_hierarchy_2__maxlevel)-1;
  }

  int find_level(Vertex_handle v) const {
    CGAL_precondition( v != Vertex_handle() );
    int level = 0;
    Vertex_handle vertex = v;
    while ( vertex->up() != Vertex_handle() ) {
      vertex = vertex->up();
      level++;
    }

    return level;
  }

  Vertex_handle
  vertex_at_level(const Vertex_handle& v, unsigned int k) const
  {
    CGAL_precondition( k <= sdg_hierarchy_2__maxlevel );

    unsigned int level = 0;
    Vertex_handle v_at_level = v;
    while ( level < k ) {
      v_at_level = v_at_level->up();
      level++;
    }
    return v_at_level;
  }

  void print_error_message() const;
};



template<class Gt, class STag, class D_S, class LTag>
std::istream& operator>>(std::istream& is,
			 Segment_Delaunay_graph_hierarchy_2<Gt,STag,D_S,LTag>&
			 sdgh)
{
  sdgh.file_input(is);
  return is;
}

template<class Gt, class STag, class D_S, class LTag>
std::ostream& operator<<(std::ostream& os,
			 const
			 Segment_Delaunay_graph_hierarchy_2<Gt,STag,D_S,LTag>&
			 sdgh)
{
  sdgh.file_output(os);
  return os;
}
			 

} //namespace CGAL


#include <CGAL/Segment_Delaunay_graph_2/Segment_Delaunay_graph_hierarchy_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_HIERARCHY_2_H
