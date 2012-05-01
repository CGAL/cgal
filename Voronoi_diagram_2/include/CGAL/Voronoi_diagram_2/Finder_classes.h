// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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

#ifndef CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H
#define CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <map>
#include <CGAL/Triangulation_utils_2.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//-------------------------------------------------------------------
//-------------------------------------------------------------------

template<class VDA>
struct Find_next_halfedge
{
  typedef Triangulation_cw_ccw_2  CW_CCW_2;

  typedef typename VDA::Delaunay_graph::Face_handle   Delaunay_face_handle;

  void operator()(const VDA* vda, const Delaunay_face_handle& f, int i,
		  Delaunay_face_handle& fnext, int& inext) const
  {
    Delaunay_face_handle fcur = f;
    int icur = i, cw_i;
    do {
      cw_i = CW_CCW_2::cw(icur);
      fnext = fcur->neighbor(cw_i);
      inext = vda->dual().tds().mirror_index(fcur, cw_i);
      fcur = fnext;
      icur = inext;
    } while ( vda->edge_rejector()(vda->dual(), fnext, inext) );
  }
};

//-------------------------------------------------------------------

template<class VDA>
struct Find_opposite_halfedge
{
  typedef typename VDA::Delaunay_graph::Face_handle   Delaunay_face_handle;
  typedef Find_next_halfedge<VDA>                     Next_halfedge;

  void operator()(const VDA* vda, const Delaunay_face_handle& f, int i,
		  Delaunay_face_handle& fopp, int& iopp) const 
  {
    Delaunay_face_handle f1;
    int i1;
    int i_mirror = vda->dual().tds().mirror_index(f, i);

    Next_halfedge()(vda, f->neighbor(i), i_mirror, f1, i1);

    fopp = f1->neighbor(i1);
    iopp = vda->dual().tds().mirror_index(f1, i1);
  }
};

//-------------------------------------------------------------------

template<class VDA>
class Find_valid_vertex
{
 public:
  typedef typename VDA::Delaunay_graph::Face_handle  Delaunay_face_handle;
  typedef std::map<Delaunay_face_handle,bool>        Delaunay_face_map;

  Delaunay_face_handle operator()(const VDA* vda,
				  const Delaunay_face_handle& f) const
  {
    CGAL_precondition( !vda->dual().is_infinite(f) );
    Delaunay_face_map fmap;
    Delaunay_face_handle fvalid;
    find_valid_vertex(vda, f, fvalid, fmap);
    CGAL_assertion( fvalid != Delaunay_face_handle() );
    CGAL_assertion( !vda->dual().is_infinite(fvalid) );
    fmap.clear();
    return fvalid;
  }

 private:
  void find_valid_vertex(const VDA* vda, const Delaunay_face_handle& cur,
			 Delaunay_face_handle& fvalid,
			 Delaunay_face_map& fmap) const
  {
    if ( fmap.find(cur) != fmap.end() ) { return; }
    fmap[cur] = true;


    bool b[3];
    for (int i = 0; i < 3; i++) {
      b[i] = !vda->edge_rejector()(vda->dual(), cur, i);
    }

    if ( b[0] || b[1] || b[2] ) {
      if ( fvalid == Delaunay_face_handle() || cur < fvalid ) {
	if ( !vda->dual().is_infinite(cur) ) {
	  fvalid = cur;
	}
      }
    }

    for (int i = 0; i < 3; i++) {
      if ( !vda->dual().is_infinite(cur->neighbor(i)) && !b[i] ) {
	find_valid_vertex(vda, cur->neighbor(i), fvalid, fmap);
      }
    }
  }

};


//-------------------------------------------------------------------
//-------------------------------------------------------------------

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_FINDER_CLASSES_H
