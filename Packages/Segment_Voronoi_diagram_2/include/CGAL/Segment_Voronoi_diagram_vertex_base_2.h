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




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <list>
#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_storage_site_2.h>
#include <CGAL/Segment_Voronoi_diagram_simple_storage_site_2.h>

#include <CGAL/Simple_container_wrapper.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template<class Gt, class Handle, class USE_SIMPLE_STORAGE_SITE_Tag>
  struct SVDVB2_Which_storage_site;

  // use the simple storage site
  template<class Gt, class Handle>
  struct SVDVB2_Which_storage_site<Gt,Handle,Tag_false>
  {
    typedef Gt         Geom_traits;
    typedef Handle     Point_handle;
    typedef Tag_false  Storage_site_tag;

    typedef
    Segment_Voronoi_diagram_simple_storage_site_2<Geom_traits,
						  Point_handle>
    Storage_site_2;
  };

  // use the full storage site
  template<class Gt, class Handle>
  struct SVDVB2_Which_storage_site<Gt,Handle,Tag_true>
  {
    typedef Gt         Geom_traits;
    typedef Handle     Point_handle;
    typedef Tag_true   Storage_site_tag;

    typedef
    Segment_Voronoi_diagram_storage_site_2<Geom_traits,Point_handle>
    Storage_site_2;
  };


} // namespace CGALi


template < class Gt,
	   class USE_FULL_STORAGE_SITE_Tag = Tag_true,
	   class Vb = Triangulation_ds_vertex_base_2<> >
class Segment_Voronoi_diagram_vertex_base_2
  : public Vb
{
private:
  typedef typename Vb::Triangulation_data_structure  DS;

  typedef std::list<typename Gt::Point_2>     PC;

public:
  // TYPES
  //------
  typedef Gt                                  Geom_traits;
  typedef Simple_container_wrapper<PC>        Point_container;
  typedef typename Point_container::iterator  Point_handle;
  typedef Vb                                  Base;
  typedef typename Gt::Site_2                 Site_2;
  typedef USE_FULL_STORAGE_SITE_Tag           Storage_site_tag;

private:
  typedef
  CGALi::SVDVB2_Which_storage_site<Geom_traits,Point_handle,
				   Storage_site_tag>
  Which_storage_site;

public:
  typedef typename Which_storage_site::Storage_site_2  Storage_site_2;

  typedef DS           Segment_Voronoi_diagram_data_structure_2;
  
  typedef typename DS::Face_handle    Face_handle;
  typedef typename DS::Vertex_handle  Vertex_handle;


  template < typename DS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<DS2>::Other   Vb2;
    typedef Storage_site_tag                               STag;
    typedef Segment_Voronoi_diagram_vertex_base_2<Gt,STag,Vb2>  Other;
  };

  
  Segment_Voronoi_diagram_vertex_base_2 () : Vb(), ss_() {}
    
  Segment_Voronoi_diagram_vertex_base_2(const Storage_site_2& ss,
					Face_handle f)
    : Vb(f), ss_(ss)  {}

  void set_site(const Storage_site_2& ss) { ss_ = ss; }

  const Storage_site_2& storage_site() const { return ss_; }
  Site_2                site()         const { return ss_.site(); }

#if 1
  // MK::ERROR: these must be removed; one may use the storage site to
  // get access to this info...
  bool is_segment() const { return ss_.is_segment(); }
  bool is_point()   const { return ss_.is_point(); }
#endif

  //the following trivial is_valid to allow
  // the user of derived face base classes 
  // to add their own purpose checking
  bool is_valid(bool /* verbose */ = false, int /* level */ = 0) const
  { return true; }

private:
  Storage_site_2 ss_;
  //  std::list<Vb>  adjseg_list; // list of adjacent segments; this is
  // important when I want to do deletions
};


CGAL_END_NAMESPACE 

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_VERTEX_BASE_2_H
