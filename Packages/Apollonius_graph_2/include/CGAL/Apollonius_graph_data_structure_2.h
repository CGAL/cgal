// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
#define CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H

#include <CGAL/Triangulation_data_structure_2.h>


CGAL_BEGIN_NAMESPACE

template <class Vb, class Fb>
class Apollonius_graph_data_structure_2
  : public Triangulation_data_structure_2<Vb, Fb>
{
public:
  typedef Apollonius_graph_data_structure_2<Vb, Fb>    Tds;
  typedef Triangulation_data_structure_2<Vb, Fb>       Tds_base;

  //  typedef typename Vb::template Rebind_TDS<Tds>::Other  Vertex_base;
  //  typedef typename Fb::template Rebind_TDS<Tds>::Other  Face_base;


public:
  typename Tds::Vertex_handle
  insert_degree_2(typename Tds::Face_handle f, int i)
  {
    /*
    // This method basically does the following transformation
    // The remove_degree_2 method performs the same operation in the
    // opposite direction
    //
    //
    //                                                *
    //                 i                             / \
    //                 *                            /   \
    //                / \                          /  f  \
    //               /   \                        / _____ \
    //              /  f  \                      / /  f1 \ \
    //             /       \                     |/   v   \|
    //  v0=ccw(i) *---------* v1=cw(i)  ===>  v0 *----*----* v1
    //             \       /                     |\   f2  /|
    //              \  g  /                      \ \_____/ /
    //               \   /                        \       /
    //                \ /                          \  g  /
    //                 *                            \   /
    //                 j                             \ /
    //                                                *
    //
    */

    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Face_handle   Face_handle;

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    Vertex_handle  v = this->create_vertex();

    Vertex_handle v0 = f->vertex( this->ccw(i) );
    Vertex_handle v1 = f->vertex( this->cw(i)  );

    Face_handle f_undef;

    Face_handle f1 = this->create_face(v0, v, v1, f_undef, f, f_undef);
    Face_handle f2 = this->create_face(v0, v1, v, f_undef, f_undef, g);

    f1->set_neighbor(0, f2);
    f1->set_neighbor(2, f2);

    f2->set_neighbor(0, f1);
    f2->set_neighbor(1, f1);

    f->set_neighbor(i, f1);
    g->set_neighbor(j, f2);

    v->set_face(f1);

    return v;
  }

  void
  remove_degree_2(typename Tds::Vertex_handle v)
  {
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Face_handle   Face_handle;

    CGAL_precondition( v->degree() == 2 );

    Face_handle f1 = v->face();
    int i = f1->index(v);

    Face_handle f2 = f1->neighbor( this->ccw(i) );
    int j = f2->index(v);

    Face_handle ff1 = f1->neighbor( i );
    Face_handle ff2 = f2->neighbor( j );

    int id1 = f1->mirror_index(i);
    int id2 = f2->mirror_index(j);

    ff1->set_neighbor(id1, ff2);
    ff2->set_neighbor(id2, ff1);

    Vertex_handle v1 = f1->vertex( this->ccw(i) );
    //    if ( v1->face() == f1 || v1->face() == f2 ) {
      v1->set_face(ff1);
      //    }

    Vertex_handle v2 = f1->vertex(  this->cw(i) );
    //    if ( v2->face() == f1 || v2->face() == f2 ) {
      v2->set_face(ff2);
      //    }

    delete_face(f1);
    delete_face(f2);

    delete_vertex(v);
  }

};




CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
