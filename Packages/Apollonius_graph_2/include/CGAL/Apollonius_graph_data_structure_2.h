// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/Apollonius_graph_data_structure_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



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
#if 0
  void
  flip(Face_handle f, int i)
  {
    CGAL_triangulation_precondition( dimension()==2 );

    // AWDG CHANGE
    // flip sould never be called if i and its mirror vertex are the same
    CGAL_triangulation_precondition( f->vertex(i) != f->mirror_vertex(i) );

    Face_handle n  = f->neighbor(i);
    // this another change: we make sure that the index of f is given
    // through the mirror vertex
    int ni = n->index( f->mirror_vertex(i) );
    // old code is here
    //  int ni = n->index(f);
    
    Vertex_handle  v_cw = f->vertex(cw(i));
    Vertex_handle  v_ccw = f->vertex(ccw(i));

    // bl == bottom left, tr == top right
    Face_handle tr = f->neighbor(ccw(i));
    Face_handle bl = n->neighbor(ccw(ni));
    int bli, tri;
    // same change here: the indices for bl and tr are given through
    // their mirror vertices
    bli = bl->index( n->mirror_vertex(ccw(ni)) );
    tri = tr->index( f->mirror_vertex(ccw(i)) );
    // old code is here
    //  bli = bl->index(n);
    //  tri = tr->index(f);
    
    f->set_vertex(cw(i), n->vertex(ni));
    n->set_vertex(cw(ni), f->vertex(i));
    
    // update the neighborhood relations
    f->set_neighbor(i, bl);
    bl->set_neighbor(bli, f);

    f->set_neighbor(ccw(i), n);
    n->set_neighbor(ccw(ni), f);
    
    n->set_neighbor(ni, tr);
    tr->set_neighbor(tri, n);
    
    if(v_cw->face() == f) {
      v_cw->set_face(n);
    }
    
    if(v_ccw->face() == n) {
      v_ccw->set_face(f);
    }
  }
#endif
  
  typename Tds::Vertex_handle
  insert_degree_2(typename Tds::Face_handle f, int i)
  {
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Face_handle   Face_handle;

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    Vertex_handle  v = create_vertex();

    Vertex_handle v0 = f->vertex(ccw(i));
    Vertex_handle v1 = f->vertex( cw(i));

    Face_handle f1 = create_face(v0, v, v1, NULL, f, NULL);
    Face_handle f2 = create_face(v0, v1, v, NULL, NULL, g);

    f1->set_neighbor(0, f2);
    f1->set_neighbor(2, f2);

    f2->set_neighbor(0, f1);
    f2->set_neighbor(1, f1);

    f->set_neighbor(i, f1);
    g->set_neighbor(j, f2);

    v->set_face(f1);

    return v;
  }

#if 0
  Vertex_handle
  insert_in_face(Face_handle f)
  {
    // AWDG CHANGE
    // new code for inserting a vertex in a face
    // the old one relies on getting right the index of a face which can
    // give a wrong answer in our case, because the index of a face is
    // not longer unique; this is because two faces can have two common edges

    CGAL_triangulation_precondition( f != NULL && dimension()== 2);
    Vertex_handle  v = create_vertex();

    Vertex_handle v0 = f->vertex(0);
    Vertex_handle v2 = f->vertex(2);
    Vertex_handle v1 = f->vertex(1);
    
    Face_handle n1 = f->neighbor(1);
    Face_handle n2 = f->neighbor(2);

    // the computation for i1 and i2 needs to be changed
    int i1(0), i2(0);
    if ( n1 != NULL ) {
      i1 = n1->index( f->mirror_vertex(1) );
    }
    if ( n2 != NULL ) {
      i2 = n2->index( f->mirror_vertex(2) );
    }

    Face_handle f1 = create_face(v0, v, v2, f, n1, NULL);
    Face_handle f2 = create_face(v0, v1, v, f, NULL, n2);

    f1->set_neighbor(2, f2);
    f2->set_neighbor(1, f1);
    if (n1 != NULL) {
      n1->set_neighbor(i1,f1);
    }
    if (n2 != NULL) {
      n2->set_neighbor(i2,f2);
    }

    f->set_vertex(0, v);
    f->set_neighbor(1, f1);
    f->set_neighbor(2, f2);

    if( v0->face() == f  ) {  v0->set_face(f2); }
    v->set_face(f);

    return v;
 }
#endif

  void
  remove_degree_2(typename Tds::Vertex_handle v)
  {
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Face_handle   Face_handle;

    CGAL_precondition( v->degree() == 2 );

    Face_handle f1 = v->face();
    int i = f1->index(v);

    Face_handle f2 = f1->neighbor( ccw(i) );
    int j = f2->index(v);

    Face_handle ff1 = f1->neighbor( i );
    Face_handle ff2 = f2->neighbor( j );

    int id1 = f1->mirror_index(i);
    int id2 = f2->mirror_index(j);

    ff1->set_neighbor(id1, ff2);
    ff2->set_neighbor(id2, ff1);

    Vertex_handle v1 = f1->vertex( ccw(i) );
    //    if ( v1->face() == f1 || v1->face() == f2 ) {
      v1->set_face(ff1);
      //    }

    Vertex_handle v2 = f1->vertex(  cw(i) );
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

    CGAL_triangulation_precondition( is_valid() );

    std::cout << "inside join: just started..." << std::endl;

    Vertex* v1 = f->vertex( ccw(i) );
    Vertex* v2 = f->vertex(  cw(i) );

    CGAL_triangulation_precondition( v == v1 || v == v2 );

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

    Face* tl = f->neighbor(  cw(i) );
    Face* tr = f->neighbor( ccw(i) );

    int itl = f->mirror_index(  cw(i) );
    int itr = f->mirror_index( ccw(i) );

    Face* bl = g->neighbor( ccw(j) );
    Face* br = g->neighbor(  cw(j) );

    int ibl = g->mirror_index( ccw(j) );
    int ibr = g->mirror_index(  cw(j) );

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

    CGAL_triangulation_postcondition( is_valid() );

    return v1;
  }

  inline
  Vertex*
  join_vertices(Face* f, int i)
  {
    return join_vertices(f, i, f->vertex( ccw(i) ));
  }
#endif
};




CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_DATA_STRUCTURE_2_H
