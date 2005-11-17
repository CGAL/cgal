// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef BSO_INTERSECTION_FUNCTOR
#define Bso_INTERSECTION_FUNCTOR

#include <CGAL/Boolean_set_operations_2/Bso_dcel.h>
#include <CGAL/Boolean_set_operations_2/Bso_base_functor.h>

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Bso_intersection_functor : public Bso_base_functor<Traits_>
{
public:

  typedef Traits_                        Traits;
  typedef Bso_dcel<Traits>              Dcel;
  typedef Arrangement_2<Traits, Dcel>    Bso_arrangement;
  typedef Bso_base_functor<Traits_>     Base; 


  typedef typename Bso_arrangement::Face_const_handle       Face_const_handle;
  typedef typename Bso_arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Bso_arrangement::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Bso_arrangement::Face_handle        Face_handle;
  typedef typename Bso_arrangement::Halfedge_handle    Halfedge_handle;
  typedef typename Bso_arrangement::Vertex_handle      Vertex_handle;



  // default constructor
  Bso_intersection_functor(Traits* tr) : Base(tr)
  {}

  void create_edge(Halfedge_const_handle h1,
                   Face_const_handle f2,
                   Halfedge_handle res_h)
  {
    if(f2->contained())
      this ->insert_edge(res_h, h1->face()->contained());
  }

  void create_edge(Face_const_handle f1,
                   Halfedge_const_handle h2,
                   Halfedge_handle res_h)
  {
    if(f1->contained())
      this ->insert_edge(res_h, h2->face()->contained());
  }

   void create_edge(Halfedge_const_handle h1,
                    Halfedge_const_handle h2,
                    Halfedge_handle res_h)
  {
    if(h1->face()->contained() && h2->face()->contained() ||
       h1->twin()->face()->contained() && h2->twin()->face()->contained())
       // in that case, the face incident to the overlap Halfedge belogns to
       //intersection face (at the result arr)

       //        _________
       //        |       |
       //        |   |---|
       //        |   |   |
       //        |   |___|
       //        |       |
       //        ---------

       this->insert_edge(res_h, h1->face()->contained());
       
  }

};


CGAL_END_NAMESPACE

#endif
