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




#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_DATA_STRUCTURE_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_DATA_STRUCTURE_2_H

#include <list>
#include <boost/tuple/tuple.hpp>
#include <CGAL/utility.h>
#include <CGAL/Apollonius_graph_data_structure_2.h>

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

CGAL_BEGIN_NAMESPACE

template <class Vb, class Fb>
class Segment_Voronoi_diagram_data_structure_2
  : public Apollonius_graph_data_structure_2<Vb, Fb>
{
private:
  typedef Apollonius_graph_data_structure_2<Vb,Fb>   Base;

public:
  typedef Segment_Voronoi_diagram_data_structure_2<Vb, Fb>    Tds;
  typedef typename Base::Tds_base                             Tds_base;

  typedef typename Tds_base::Vertex_handle    Vertex_handle;
  typedef typename Tds_base::Face_handle      Face_handle;
  typedef typename Tds_base::Face_circulator  Face_circulator;

  typedef typename Vb::template Rebind_TDS<Tds>::Other  Vertex_base;
  typedef typename Fb::template Rebind_TDS<Tds>::Other  Face_base;


  /*
  // The following method preforms a split operation of the vertex v
  // using the faces f1 and g1. The split operation is shown
  // below.
  // The names of the variables in the method correspond to the
  // quantities in the drawings below
  //
  // The configuration before the split:
  //
  //                  cw(i1)   v3   ccw(i2)
  //                     *-----*-----*
  //                    / \    |    / \
  //                   /   \ f1|f2 /   \
  //                  /     \  |  /     \
  //                 /       \ | /       \
  //                /         \|/v        \
  //               *-----------*-----------*
  //                \         /|\         /
  //                 \       / | \       /
  //                  \     /  |  \     /
  //                   \   / g2|g1 \   /
  //                    \ /    |    \ /
  //                     *-----*-----*
  //                 ccw(j2)   v4   cw(j1)
  //
  //
  // The configuration after the split:
  //
  //
  //               cw(i1)      v3     ccw(i2)
  //                 *---------*---------*
  //                / \       / \       / \
  //               /   \  f1 /   \  f2 /   \
  //              /     \   /  f  \   /     \
  //             /       \ /     v2\ /       \
  //            *---------*---------*---------*
  //             \       / \v1     / \       /
  //              \     /   \  g  /   \     /
  //               \   /  g2 \   /  g1 \   /
  //                \ /       \ /       \ /
  //                 *---------*---------*
  //              ccw(j2)      v4      cw(j1)
  //
  */
//  Quadruple<Vertex_handle, Vertex_handle, Face_handle, Face_handle>
  boost::tuples::tuple<Vertex_handle,Vertex_handle,Face_handle,Face_handle>
  split_vertex(Vertex_handle v, Face_handle f1, Face_handle g1)
  {
    CGAL_precondition( this->dimension() == 2 );
    CGAL_precondition( f1 != Face_handle() && f1->has_vertex(v) );
    CGAL_precondition( g1 != Face_handle() && g1->has_vertex(v) );

    typedef Vertex_handle VH;
    typedef Face_handle   FH;

    // 1. first we read some information that we will need
    int i1 = f1->index(v);
    int j1 = g1->index(v);
    Face_handle f2 = f1->neighbor( this->cw(i1) );
    Face_handle g2 = g1->neighbor( this->cw(j1) );

    int i2 = f2->index(v);
    int j2 = g2->index(v);

    Vertex_handle v3 = f1->vertex( this->ccw(i1) );
    Vertex_handle v4 = g1->vertex( this->ccw(j1) );

    // lst is the list of faces adjecent to v stored in
    // counterclockwise order from g2 to f1) inclusive.
    // the list idx contains the indices of v in the
    // faces in lst.
    std::list<Face_handle> lst;
    std::list<int>         idx;

    Face_circulator fc(v, g1);
    Face_handle ff(fc);
    while ( ff != f2 ) {
      lst.push_back( ff );
      idx.push_back( ff->index(v) );
      fc++;
      ff = Face_handle(fc);
    }
    lst.push_back( ff );
    idx.push_back( ff->index(v) );

    // 2. we create the new vertices and the two new faces
    Vertex_handle v1 = v;
    Vertex_handle v2 = this->create_vertex();
    Face_handle f = this->create_face(v1, v2, v3);
    Face_handle g = this->create_face(v2, v1, v4);

    // 3. we update the adjacency information for the new vertices and
    //    the new faces
    f->set_neighbor(0, f2);
    f->set_neighbor(1, f1);
    f->set_neighbor(2, g);
    g->set_neighbor(0, g2);
    g->set_neighbor(1, g1);
    g->set_neighbor(2, g);
    v1->set_face(f);
    v2->set_face(g);

    // 4. update the vertex for the faces f2 through g1 in
    //    counterclockwise order
    typename std::list<Face_handle>::iterator fit = lst.begin();
    typename std::list<int>::iterator         iit = idx.begin();
    for (; fit != lst.end(); ++fit, ++iit) {
      (*fit)->set_vertex(*iit, v2);
    }

    lst.clear();
    idx.clear();

    // 5. make f and g the new neighbors of f1, f2 and g1, g2
    //    respectively.
    f1->set_neighbor( this->cw(i1), f );
    f2->set_neighbor(this->ccw(i2), f );
    g1->set_neighbor( this->cw(j1), g );
    g2->set_neighbor(this->ccw(j2), g );

    // 6. return the new stuff
    return boost::tuples::tuple<VH, VH, FH, FH>(v1, v2, f, g);
//    return Quadruple<VH, VH, FH, FH>(v1, v2, f, g);
  }
};


CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_DATA_STRUCTURE_2_H
