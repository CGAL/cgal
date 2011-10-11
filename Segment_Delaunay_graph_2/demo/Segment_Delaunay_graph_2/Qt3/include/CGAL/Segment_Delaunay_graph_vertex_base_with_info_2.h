// Copyright (c) 2003,2004,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_WITH_INFO_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_WITH_INFO_2_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>

#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_site_2.h>
#include <CGAL/Segment_Delaunay_graph_simple_storage_site_2.h>

namespace CGAL {


template < class Vbb, class Info_ >
class Segment_Delaunay_graph_vertex_base_with_info_2
  : public Vbb
{
public:
  // TYPES
  //------
  typedef Vbb                                      Base;
  typedef Info_                                    Info;

  typedef typename Base::Storage_traits            Storage_traits;
  typedef typename Base::Geom_traits               Geom_traits;
  typedef typename Base::Site_2                    Site_2;
  typedef typename Base::Storage_site_2            Storage_site_2;
  typedef typename Base::Data_structure            Data_structure;

  typedef typename Data_structure::Face_handle     Face_handle;
  typedef typename Data_structure::Vertex_handle   Vertex_handle;


  template < typename DS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<DS2>::Other             Vb2;
    typedef Segment_Delaunay_graph_vertex_base_with_info_2<Vb2,Info> Other;
  };


  Segment_Delaunay_graph_vertex_base_with_info_2 () : Vbb(), info_() {}

  Segment_Delaunay_graph_vertex_base_with_info_2(const Storage_site_2& ss,
						 Face_handle f)
    : Vbb(ss, f), info_()  {}

  void set_info(const Info& info) { info_ = info; }
  const Info& info() const { return info_; }

private:
  Info info_;
};


} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_VERTEX_BASE_WITH_INFO_2_H
