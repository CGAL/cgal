// Copyright (c) 2003,2004,2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_site_2.h>
#include <CGAL/Segment_Delaunay_graph_simple_storage_site_2.h>



namespace CGAL {

template < class STraits, class Vb = Triangulation_ds_vertex_base_2<> >
class Segment_Delaunay_graph_vertex_base_2
  : public Vb
{
private:
  typedef typename Vb::Triangulation_data_structure  D_S;
  typedef Vb                                         Base;

public:
  // TYPES
  //------
  typedef STraits                                  Storage_traits;
  typedef typename Storage_traits::Geom_traits     Geom_traits;
  typedef typename Geom_traits::Point_2            Point;
  typedef typename Geom_traits::Site_2             Site_2;
  typedef typename Storage_traits::Storage_site_2  Storage_site_2;
  typedef D_S                                       Data_structure;
  
  typedef typename D_S::Face_handle                 Face_handle;
  typedef typename D_S::Vertex_handle               Vertex_handle;


  template < typename DS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<DS2>::Other       Vb2;
    typedef Segment_Delaunay_graph_vertex_base_2<STraits,Vb2>  Other;
  };

  
  Segment_Delaunay_graph_vertex_base_2 () : Vb(), ss_() {}
    
  Segment_Delaunay_graph_vertex_base_2(const Storage_site_2& ss,
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


} //namespace CGAL 

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_2_H
