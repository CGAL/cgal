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

#if 0
  inline
  Vertex*
  join_vertices(Face* f, int i, Vertex* v)
  {
    // this methods does the "join"-operation and preserves
    // the vertex v among the two vertices that define the edge (f, i) 

    CGAL_precondition( is_valid() );

    std::cout << "inside join: just started..." << std::endl;

    Vertex* v1 = f->vertex( this->ccw(i) );
    Vertex* v2 = f->vertex( this->cw(i)  );

    CGAL_precondition( v == v1 || v == v2 );

    if ( v == v2 ) {
      return join_vertices(f->neighbor(i), f->mirror_index(i), v);
    }

    int deg2 = v2->degree();

    std::cout << "degrees 1, 2: " << v1->degree()
	      << " " << deg2 << std::endl;

    if ( deg2 == 3 ) {
      remove_degree_3(v2, f);
      return v1;
    } else if ( deg2 == 2 ) {
      remove_degree_2(v2);
      return v1;
    }

    /*
    // The following drawing corrsponds to the variables
    // used in this part...
    // The vertex v1 is returned...
    //
    //      itl       i=v3      itr
    //       *---------*---------*
    //        \       / \       /
    //         \  tl /   \  tr /
    //          \   /  f  \   /
    //           \ /       \ /
    //  v1=ccw(i) *---------*  cw(i)=v2
    //           / \       / \
    //          /   \  g  /   \
    //         /  bl \   /  br \
    //        /       \ /       \
    //       *---------*---------*
    //      ibl       j=v4      ibr
    //                                                           
    // The situation after the "join"-operation is as follows:
    //
    //                 i
    //           *-----*-----*
    //            \    |    /
    //             \ tl|tr /
    //              \  |  /
    //               \ | /
    //                \|/
    //                 *  v1
    //                /|\
    //               / | \
    //              /  |  \
    //             / bl|br \
    //            /    |    \
    //           *-----*-----*
    //
    */

    std::cout << "inside join: degree > 3" << std::endl;

    // first we register all the needed info
    Face* g = f->neighbor(i);
    int j = f->mirror_index(i);

    Face* tl = f->neighbor( this->cw(i)  );
    Face* tr = f->neighbor( this->ccw(i) );

    int itl = f->mirror_index( this->cw(i)  );
    int itr = f->mirror_index( this->ccw(i) );

    Face* bl = g->neighbor( this->ccw(j) );
    Face* br = g->neighbor( this->cw(j)  );

    int ibl = g->mirror_index( this->ccw(j) );
    int ibr = g->mirror_index( this->cw(j)  );

    // we need to store the faces adjacent to v2 as well as the
    // indices of v2 w.r.t. these faces, so that afterwards we can set 
    // v1 to be the vertex for these faces
    std::vector<Face*> star_faces_of_v2;
    std::vector<int> star_indices_of_v2;
    typename Tds_base::Face_circulator fc_start(v2);
    typename Tds_base::Face_circulator fc = fc_start;

    std::cout << "inside join..." << std::endl;

    do {
      Face* ff = &(*fc);
      star_faces_of_v2.push_back(ff);
      star_indices_of_v2.push_back(ff->index(v2));
      ++fc;
    } while ( fc != fc_start );

    CGAL_assertion( int(star_faces_of_v2.size()) == deg2 );

    // from this point and on we modify the values

    std::cout << "inside join: changing info..." << std::endl;

    // first set the neighbors
    tl->set_neighbor(itl, tr);
    tr->set_neighbor(itr, tl);

    bl->set_neighbor(ibl, br);
    br->set_neighbor(ibr, bl);

    // make sure that all the faces containing v2 as a vertex, now
    // contain v1
    for (unsigned int k = 0; k < star_faces_of_v2.size(); k++) {
      int id = star_indices_of_v2[k];
      CGAL_assertion( star_faces_of_v2[k]->vertex(id) == v2 );
      star_faces_of_v2[k]->set_vertex( id, v1 );
    }

    // then make sure that all the vertices have correct pointers to 
    // faces
    Vertex* v3 = f->vertex(i);
    Vertex* v4 = g->vertex(j);
    if ( v3->face() == f )  v3->set_face(tr);
    if ( v4->face() == g )  v4->set_face(br);
    if ( v1->face() == f || v1->face() == g ) v1->set_face(tl);


    typename Tds_base::Face_iterator fit = faces_begin();
    for (; fit != faces_end(); ++fit) {
      int id;
      CGAL_assertion( !fit->has_vertex(v2, id) );
    }

    // memory management
    star_faces_of_v2.clear();
    star_indices_of_v2.clear();

    delete_face(f);
    delete_face(g);

    delete_vertex(v2);

    CGAL_postcondition( is_valid() );

    return v1;
  }

  inline
  Vertex*
  join_vertices(Face* f, int i)
  {
    return join_vertices(f, i, f->vertex( this->ccw(i) ));
  }
#endif
};




CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
