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
#include <CGAL/Segment_Voronoi_diagram_storage_site_2.h>

CGAL_BEGIN_NAMESPACE


template < class Gt, class PointHandle,
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
  typedef PointHandle             Point_handle;
  typedef Vb                      Base;
  typedef typename Gt::Site_2     Site_2;

#ifdef USE_STORAGE_SITE
  typedef
  Segment_Voronoi_diagram_storage_site_2<Gt,Point_handle>
  Storage_site_2;
#endif

  typedef SVDDS           Segment_Voronoi_diagram_data_structure_2;
  
  typedef typename SVDDS::Face_handle    Face_handle;
  typedef typename SVDDS::Vertex_handle  Vertex_handle;


  template < typename SVDDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<SVDDS2>::Other  Vb2;
    //    typedef Segment_Voronoi_diagram_vertex_base_2<Site_2,Vb2> Other;
    typedef
    Segment_Voronoi_diagram_vertex_base_2<Gt,Point_handle,Vb2>  Other;
  };

  
#ifdef USE_STORAGE_SITE
  Segment_Voronoi_diagram_vertex_base_2 () : Vb(), ss_() {}
    
  Segment_Voronoi_diagram_vertex_base_2(const Storage_site_2& ss,
					Face_handle f)
    : Vb(f), ss_(ss)  {}
#else
  Segment_Voronoi_diagram_vertex_base_2 () : Vb(), s_() {}
    
  Segment_Voronoi_diagram_vertex_base_2(const Site_2& s,
					Face_handle f)
    : Vb(f), s_(s)  {}
#endif

#ifdef USE_STORAGE_SITE
  void set_site(const Storage_site_2& ss) {
    ss_ = ss;
  }
#else
  void set_site(const Site_2& s) {
    s_ = s;
  }  
#endif

#ifdef USE_STORAGE_SITE
  const Storage_site_2& storage_site() const { return ss_; }
  Site_2                site()         const { return ss_.site(); }

  bool is_segment() const { return ss_.is_segment(); }
  bool is_point()   const { return ss_.is_point(); }
#else
  const Site_2&         site()         const { return s_; }

  bool is_segment() const { return s_.is_segment(); }
  bool is_point()   const { return s_.is_point(); }
#endif

  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
  { return true; }

private:
#ifdef USE_STORAGE_SITE
  Storage_site_2 ss_;
#else
  Site_2 s_;
#endif
  //  std::list<Vb>  adjseg_list; // list of adjacent segments; this is
  // important when I want to do deletions
  //  bool _is_point; // false if it is a segment; true otherwise
  // to be used for intersecting segments:
  //  bool is_simple;
};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H
