// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Segment_Voronoi_diagram_vertex_base_2.h
// package       : Segment_Voronoi_diagram_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H


#include <CGAL/Triangulation_ds_vertex_base_2.h>

CGAL_BEGIN_NAMESPACE


template < class Gt,
	   class Vb = Triangulation_ds_vertex_base_2<> >
class Segment_Voronoi_diagram_vertex_base_2
  : public Vb
{
private:
  typedef typename Vb::Triangulation_data_structure    SVDDS;
public:
  // TYPES
  //------
  typedef Gt                      Geom_traits;
  typedef Vb                      Base;
  typedef typename Gt::Site_2     Site_2;
  typedef typename Gt::Point_2    Point_2;
  typedef typename Gt::Segment_2  Segment_2;

  typedef SVDDS           Segment_Voronoi_diagram_data_structure_2;
  
  typedef typename SVDDS::Face_handle    Face_handle;
  typedef typename SVDDS::Vertex_handle  Vertex_handle;


  template < typename SVDDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<SVDDS2>::Other  Vb2;
    typedef Segment_Voronoi_diagram_vertex_base_2<Gt,Vb2>    Other;
  };

  
  Segment_Voronoi_diagram_vertex_base_2 () : Vb(), _s() {}
    
#if 0
  Segment_Voronoi_diagram_vertex_base_2(const Point_2 & p,
					Face_handle f = Face_handle(NULL)) 
    :  _s(p), _f(f), is_point(false)  {}

  Segment_Voronoi_diagram_vertex_base_2(const Segment_2 & s,
					Face_handle f = Face_handle(NULL))
    :  _s(s), _f(f), is_point(false)  {}
#endif

  Segment_Voronoi_diagram_vertex_base_2(const Site_2 & t, Face_handle f)
    :  Vb(f), _s(t)  {}


  inline void set_point(const Point_2& p)
    { _s.set_point(p); }
  inline void set_segment(const Segment_2& s)
    { _s.set_segment(s); }
  inline void set_site(const Site_2& t)
    { _s = t; }
  
  //  inline void set_face(Face_handle f) { _f = f; }
 
  Point_2    point() const
  { CGAL_precondition( is_point() ); return _s.point(); }

  Segment_2  segment() const
  { CGAL_precondition( is_segment() ); return _s.segment(); }

  const Site_2& site() const { return _s; }

  Point_2   source() const
  { CGAL_precondition( is_segment() ); return _s.source(); }

  Point_2   target() const
  { CGAL_precondition( is_segment() ); return _s.target(); }
 
#if 0
  inline Point_2&    point()
    { CGAL_precondition( is_point() ); return _s.point(); }
  inline Segment_2&  segment()
    { CGAL_precondition( is_segment() ); return _s.segment(); }
  inline Site_2&     site() { return _s; }
#endif

  bool is_segment() const { return _s.is_segment(); }
  bool is_point()   const { return _s.is_point(); }

   
  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
  { return true; }

private:
  Site_2 _s;
  //  std::list<Vb>  adjseg_list; // list of adjacent segments; this is
  // important when I want to do deletions
  //  bool _is_point; // false if it is a segment; true otherwise
  // to be used for intersecting segments:
  //  bool is_simple;
};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H
