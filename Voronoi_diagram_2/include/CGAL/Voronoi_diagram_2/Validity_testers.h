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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_VALIDITY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_VALIDITY_TESTERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <algorithm>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Voronoi_diagram_2/Finder_classes.h>


namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================


template<class VDA, class Base_it>
class Edge_validity_tester
{
  // tests whether a halfedge has as face a face with zero area.
 private:
  const VDA* vda_;

 private:
  typedef Triangulation_cw_ccw_2                       CW_CCW_2;
  // Base_it is essentially VDA::Edges_iterator_base
  typedef Base_it                                      Edges_iterator_base;
  typedef typename VDA::Halfedge_handle                Halfedge_handle;
  typedef typename VDA::Delaunay_graph::Vertex_handle  Delaunay_vertex_handle;

 public:
  Edge_validity_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Edges_iterator_base& eit) const {
    CGAL_assertion( !vda_->edge_rejector()(vda_->dual(), eit->dual()) );

    int cw_i = CW_CCW_2::cw( eit->dual().second );
    CGAL_assertion_code( int ccw_i = CW_CCW_2::ccw( eit->dual().second ); )

    CGAL_assertion_code(Delaunay_vertex_handle v_ccw_i = eit->dual().first->vertex(ccw_i);)
    CGAL_assertion(  !vda_->face_rejector()(vda_->dual(), v_ccw_i)  );

    Delaunay_vertex_handle v_cw_i = eit->dual().first->vertex(cw_i);
    if (  !vda_->face_rejector()(vda_->dual(), v_cw_i)  ) {
      return false;
    }

    Halfedge_handle he(eit);
    Halfedge_handle he_opp = eit->opposite();

    CGAL_assertion( he_opp->opposite() == he );

    return he->face()->dual() < he_opp->face()->dual();
  }
};


//=========================================================================
//=========================================================================


template<class VDA>
class Vertex_validity_tester
{
 private:
  const VDA* vda_;

 private:
  typedef typename VDA::Delaunay_graph::Face_handle     Delaunay_face_handle;
  typedef typename VDA::Delaunay_graph::Finite_faces_iterator
  Delaunay_faces_iterator;

 public:
  Vertex_validity_tester(const VDA* vda = NULL) : vda_(vda) {}

  bool operator()(const Delaunay_faces_iterator& fit) const {
    Delaunay_face_handle f(fit);
    Delaunay_face_handle fvalid = Find_valid_vertex<VDA>()(vda_,f);
    return f != fvalid;
  }
};


//=========================================================================
//=========================================================================

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_VALIDITY_TESTERS_H
