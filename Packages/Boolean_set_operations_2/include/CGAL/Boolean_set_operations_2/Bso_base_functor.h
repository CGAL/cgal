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

#ifndef BSO_BASE_FUNCTOR
#define BSO_BASE_FUNCTOR

#include <CGAL/Boolean_set_operations_2/Bso_dcel.h>
#include <CGAL/Unique_hash_map.h> 
#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_accessor.h>

CGAL_BEGIN_NAMESPACE


template <class Traits_>
class Bso_base_functor
{
public:

  typedef Traits_                           Traits;
  typedef Bso_dcel<Traits>                 Dcel;
  typedef Arrangement_2<Traits, Dcel>       Bso_arrangement;

  typedef typename Bso_arrangement::Face_const_handle       Face_const_handle;
  typedef typename Bso_arrangement::Vertex_const_handle     Vertex_const_handle;
  typedef typename Bso_arrangement::Halfedge_const_handle   Halfedge_const_handle;

  typedef typename Bso_arrangement::Face_handle        Face_handle;
  typedef typename Bso_arrangement::Halfedge_handle    Halfedge_handle;
  typedef typename Bso_arrangement::Vertex_handle      Vertex_handle;
  
  typedef Unique_hash_map<Vertex_handle, Vertex_handle> Vertex_map;



  // default constructor
  Bso_base_functor(Traits *tr) : m_vertices_map(Vertex_handle()),
                                  m_res_arr(new Bso_arrangement(tr)),
                                  m_arr_access (*m_res_arr)
  {}

   void create_face (Face_const_handle f1,
                     Face_const_handle f2,
                     Face_handle res_f)
  {}

  void create_vertex(Halfedge_const_handle h1,
                     Halfedge_const_handle h2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Vertex_const_handle v2,
                     Vertex_handle  res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Halfedge_const_handle h2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Halfedge_const_handle h1,
                     Vertex_const_handle v2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Face_const_handle f1,
                     Vertex_const_handle v2,
                     Vertex_handle res_v)
  {}

  void create_vertex(Vertex_const_handle v1,
                     Face_const_handle f2,
                     Vertex_handle res_v)
  {}

  void create_edge(Halfedge_const_handle h1,
                   Halfedge_const_handle h2,
                   Halfedge_handle res_h)
  {}

  void create_edge(Halfedge_const_handle h1,
                   Face_const_handle f2,
                   Halfedge_handle res_h)
  {}

  void create_edge(Face_const_handle f1,
                   Halfedge_const_handle h2,
                   Halfedge_handle res_h)
  {}


  void insert_edge(Halfedge_handle he, bool is_face_contained)
  {
    Vertex_handle s = he->source();
    Vertex_handle t = he->target();

    Vertex_handle s_in = m_vertices_map[s];
    Vertex_handle t_in = m_vertices_map[t];

    Vertex_handle def_v;
    bool s_exist = (s_in != def_v);
    bool t_exist = (t_in != def_v);

    CGAL_precondition_code(typename Traits::Compare_xy_2  cmp_xy =
      m_res_arr->get_traits()->compare_xy_2_object(););
    CGAL_assertion(cmp_xy(s->point(), t->point()) == LARGER);


    if(!s_exist && !t_exist)
    {
      //insert in face interior (allways the unbounded face)
      Halfedge_handle res = 
        m_res_arr->insert_in_face_interior(he->curve(),
                                          m_res_arr->unbounded_face());
      
      // res is directed from left to right (opposite to 'he')
      Vertex_handle s2 = res->source();
      Vertex_handle t2 = res->target();

      m_vertices_map[s] = t2;
      m_vertices_map[t] = s2;
      
    }
    else 
      if(s_exist && !t_exist)
      {
        Halfedge_handle res = 
          m_res_arr->insert_from_right_vertex(he->curve(), s_in);
      
        // res is directed from right to left
        Vertex_handle t2 = res->target();

        m_vertices_map[t] = t2;
       
      }
    else 
      if(!s_exist && t_exist)
      {
        Halfedge_handle res =
          m_res_arr->insert_from_left_vertex(he->curve(), t_in);

        //res is firected from left to right
        

        Vertex_handle t2 = res->target();

        m_vertices_map[s] = t2;
      }
      else
      {
        CGAL_assertion(s_exist && t_exist);
        Halfedge_handle res;
        if(s_in->degree() == 1 && t_in->degree() == 1)
        {
          bool new_face_created;
          res = m_arr_access.insert_at_vertices_ex (he->curve(),
                                                    s_in->incident_halfedges(),
                                                    t_in->incident_halfedges(),
                                                    LARGER,
                                                    new_face_created);
          
          if (new_face_created)
          {
            // In case a new face has been created (pointed by the new halfedge
            // we obtained), we have to examine the holes and isolated vertices
            // in the existing face (pointed be the twin halfedge) and relocate
            // the relevant features in the new face.
            m_arr_access.relocate_in_new_face (res);
          }
        }
        else
        {
          res = m_res_arr->insert_at_vertices(he->curve(), s_in, t_in);
        }

       
        //res = m_res_arr->insert_at_vertices(he->curve(), s_in, t_in);
        CGAL_assertion(res->direction() == LARGER);
        if(res->face() != res->twin()->face()) // new face created (update it)
        {
          //res must be incident to the new face
          res->face()->set_contained(is_face_contained);
        }
      }
    }

 
  Bso_arrangement* result_arr()
  {
    return m_res_arr;
  }


  

protected:

  Vertex_map                      m_vertices_map;
  Bso_arrangement*               m_res_arr;
  Arr_accessor<Bso_arrangement>  m_arr_access;

};


CGAL_END_NAMESPACE

#endif
