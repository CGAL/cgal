// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_OVERLAY_FUNCTOR_H
#define CGAL_ENVELOPE_OVERLAY_FUNCTOR_H

#include <CGAL/license/Envelope_3.h>


#include <iostream>

namespace CGAL {

template < class MinimizationDiagram_2>
class Envelope_overlay_functor
{
public:
  typedef MinimizationDiagram_2                                  Minimization_diagram_2;

  typedef typename Minimization_diagram_2::Face_const_handle     Face_const_handle1;
  typedef typename Minimization_diagram_2::Face_const_handle     Face_const_handle2;

  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_const_handle1;
  typedef typename Minimization_diagram_2::Vertex_const_handle   Vertex_const_handle2;

  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_const_handle1;
  typedef typename Minimization_diagram_2::Halfedge_const_handle Halfedge_const_handle2;

  typedef typename Minimization_diagram_2::Face_handle           Face_handle;
  typedef typename Minimization_diagram_2::Vertex_handle         Vertex_handle;
  typedef typename Minimization_diagram_2::Halfedge_handle       Halfedge_handle;

  typedef typename Minimization_diagram_2::Face_handle           Res_face_handle;
  typedef typename Minimization_diagram_2::Halfedge_handle       Res_halfedge_handle;
  typedef typename Minimization_diagram_2::Vertex_handle         Res_vertex_handle;

protected:

  typedef typename Minimization_diagram_2::Dcel::Dcel_data_iterator Envelope_data_iterator;
  typedef std::pair<Vertex_handle, Face_handle>                  Vertex_face_pair;
  struct Less_vertex_face_pair
  {
    bool operator() (const Vertex_face_pair& vf1,
                     const Vertex_face_pair& vf2) const
    {
      Vertex_handle v1 = vf1.first, v2 = vf2.first;
      Face_handle f1 = vf1.second, f2= vf2.second;
      return (&*v1 < &*v2 ||
              (&*v1 == &*v2 && &*f1 < &*f2));
    }
  };
  typedef std::map<Vertex_face_pair, Halfedge_handle,
                   Less_vertex_face_pair>                        Boundary_cache;

public:

  Envelope_overlay_functor(Minimization_diagram_2& md1,
                           Minimization_diagram_2& md2,
                           Minimization_diagram_2& result)
    : m_1(md1), m_2(md2), m_result(result)
  {
  }

  ~Envelope_overlay_functor()
  {
    traversed_vertices.clear();
  }

  void create_face (Face_const_handle1 f1, Face_const_handle2 f2, Res_face_handle res_f)
  {
    res_f->set_aux_source(0, m_1.non_const_handle(f1));
    res_f->set_aux_source(1, m_2.non_const_handle(f2));
  }

  void create_vertex(Halfedge_const_handle1 h1,
                     Halfedge_const_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(h1));
    res_v->set_aux_source(1, m_2.non_const_handle(h2));
    //res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    //res_v->set_is_intersection(false);

    if (v1->is_isolated() && v2->is_isolated())
    {
      res_v->set_is_equal_aux_data_in_face(0, v1->get_is_equal_data_in_face());
      res_v->set_is_equal_aux_data_in_face(1, v2->get_is_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(0, v1->get_has_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(1, v2->get_has_equal_data_in_face());
    }
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Halfedge_const_handle2 h2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(h2));
    //res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Halfedge_const_handle1 h1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(h1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    //res_v->set_is_intersection(true);
    // res_v cannot be isolated
  }

  void create_vertex(Face_const_handle1 f1,
                     Vertex_const_handle2 v2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(f1));
    res_v->set_aux_source(1, m_2.non_const_handle(v2));
    //res_v->set_is_intersection(false);

    if (v2->is_isolated())
    {
      // the res_v is also isolated, and we should update the is_equal/has_equal
      // data in face information
      res_v->set_is_equal_aux_data_in_face(0, true);
      res_v->set_is_equal_aux_data_in_face(1, v2->get_is_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
      res_v->set_has_equal_aux_data_in_face(1, v2->get_has_equal_data_in_face());
    }
  }

  void create_vertex(Vertex_const_handle1 v1,
                     Face_const_handle2 f2,
                     Res_vertex_handle res_v)
  {
    res_v->set_aux_source(0, m_1.non_const_handle(v1));
    res_v->set_aux_source(1, m_2.non_const_handle(f2));
    //res_v->set_is_intersection(false);

    if (v1->is_isolated())
    {
      // the res_v is also isolated, and we should update the is_equal/has_equal
      // data in face information
      res_v->set_is_equal_aux_data_in_face(0, v1->get_is_equal_data_in_face());
      res_v->set_is_equal_aux_data_in_face(1, true);
      res_v->set_has_equal_aux_data_in_face(0, v1->get_has_equal_data_in_face());
      res_v->set_has_equal_aux_data_in_face(1, !f2->has_no_data());
    }
  }

  void create_edge(Halfedge_const_handle1 h1,
                   Halfedge_const_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    // update source
    res_h->set_aux_source(0, m_1.non_const_handle(h1));
    res_h->set_aux_source(1, m_2.non_const_handle(h2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(h1->twin()));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(h2->twin()));

    // update is_equal/has_equal data in face
    res_h->set_is_equal_aux_data_in_face(0, h1->get_is_equal_data_in_face());
    res_h->set_is_equal_aux_data_in_face(1, h2->get_is_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(0, h1->get_has_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(1, h2->get_has_equal_data_in_face());

    res_h->twin()->set_is_equal_aux_data_in_face(0, h1->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_is_equal_aux_data_in_face(1, h2->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(0, h1->twin()->get_has_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(1, h2->twin()->get_has_equal_data_in_face());

    // update is_equal/has_equal data in target
    update_halfedge_flags_on_edge(res_h, m_1.non_const_handle(h1), 0);

    // update aux_data(1)
    update_halfedge_flags_on_edge(res_h, m_2.non_const_handle(h2), 1);

    // update is_equal/has_equal data in source
    // update aux_data(0)
          update_halfedge_flags_on_edge(res_h->twin(), m_1.non_const_handle(h1->twin()), 0);

    // update aux_data(1)
          update_halfedge_flags_on_edge(res_h->twin(), m_2.non_const_handle(h2->twin()), 1);

  }

  void create_edge(Halfedge_const_handle1 h1,
                   Face_const_handle2 f2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, m_1.non_const_handle(h1));
    res_h->set_aux_source(1, m_2.non_const_handle(f2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(h1->twin()));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(f2));

    // update is_equal/has_equal data in face
    res_h->set_is_equal_aux_data_in_face(0, h1->get_is_equal_data_in_face());
    res_h->set_is_equal_aux_data_in_face(1, true);
    res_h->set_has_equal_aux_data_in_face(0, h1->get_has_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(1, !f2->has_no_data());

    res_h->twin()->set_is_equal_aux_data_in_face(0, h1->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_is_equal_aux_data_in_face(1, true);
    res_h->twin()->set_has_equal_aux_data_in_face(0, h1->twin()->get_has_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(1, !f2->has_no_data());

    // update is_equal/has_equal data in target for the first source map
          update_halfedge_flags_on_edge(res_h, m_1.non_const_handle(h1), 0);

    // update source
          update_halfedge_flags_on_edge(res_h->twin(), m_1.non_const_handle(h1->twin()), 0);

    // update is_equal/has_equal data in target for the second source map
    update_halfedge_flags_in_face(res_h, m_2.non_const_handle(f2), 1);

    // update is_equal/has_equal data in source for the second source map
    update_halfedge_flags_in_face(res_h->twin(), m_2.non_const_handle(f2), 1);
  }


  void create_edge(Face_const_handle1 f1,
                   Halfedge_const_handle2 h2,
                   Res_halfedge_handle res_h)
  {
    res_h->set_aux_source(0, m_1.non_const_handle(f1));
    res_h->set_aux_source(1, m_2.non_const_handle(h2));

    res_h->twin()->set_aux_source(0, m_1.non_const_handle(f1));
    res_h->twin()->set_aux_source(1, m_2.non_const_handle(h2->twin()));

    // update halfedge-face flags of the new halfedge
    res_h->set_is_equal_aux_data_in_face(0, true);
    res_h->set_is_equal_aux_data_in_face(1, h2->get_is_equal_data_in_face());
    res_h->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
    res_h->set_has_equal_aux_data_in_face(1, h2->get_has_equal_data_in_face());

    res_h->twin()->set_is_equal_aux_data_in_face(0, true);
    res_h->twin()->set_is_equal_aux_data_in_face(1, h2->twin()->get_is_equal_data_in_face());
    res_h->twin()->set_has_equal_aux_data_in_face(0, !f1->has_no_data());
    res_h->twin()->set_has_equal_aux_data_in_face(1, h2->twin()->get_has_equal_data_in_face());

    // update is_equal/has_equal data in target for the second source map
          update_halfedge_flags_on_edge(res_h, m_2.non_const_handle(h2), 1);

    // update source
          update_halfedge_flags_on_edge(res_h->twin(), m_2.non_const_handle(h2->twin()), 1);

    // update is_equal/has_equal data in target for the first source map
    update_halfedge_flags_in_face(res_h, m_1.non_const_handle(f1), 0);
    // update is_equal/has_equal data in source for the first source map
    update_halfedge_flags_in_face(res_h->twin(), m_1.non_const_handle(f1), 0);
  }

protected:

  template <class Halfedge_handle_t>
  void copy_halfedge_target_info(Halfedge_handle_t from, Res_halfedge_handle to, unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_target());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_target());
  }
  void set_halfedge_target_info(Res_halfedge_handle to, unsigned int id, bool info)
  {
    to->set_is_equal_aux_data_in_target(id, info);
    to->set_has_equal_aux_data_in_target(id, info);
  }

  template <class Halfedge_handle_t>
  void copy_halfedge_target_info_from_halfedge_face_info(Halfedge_handle_t from,
                                                         Res_halfedge_handle to,
                                                         unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_face());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_face());
    to->set_has_equal_aux_data_in_target_and_face(id, from->get_has_equal_data_in_face());
  }
  template <class Vertex_handle_t>
  void copy_halfedge_target_info_from_vertex_face_info(Vertex_handle_t from,
                                                       Res_halfedge_handle to,
                                                       unsigned int id)
  {
    to->set_is_equal_aux_data_in_target(id, from->get_is_equal_data_in_face());
    to->set_has_equal_aux_data_in_target(id, from->get_has_equal_data_in_face());
  }

  // find a halfedge that v is its target and f is its face
  Halfedge_handle find_halfedge_by_vertex_and_face(Vertex_handle v, Face_handle f)
  {
    // should always invoke this method when v is on the boundary of f

    // for the complexity of the total algorithm, we only loop over
    // the halfedges of each vertex once and cache the triples of
    // vertex-face-halfedge for future such questions
    Vertex_face_pair query(v, f);
    typename Boundary_cache::iterator iter = traversed_vertices.find(query);
    Halfedge_handle result;
    if (iter == traversed_vertices.end())
    {
      // first time to check this vertex - traverse all its halfedges
      // and update the map
      typename Minimization_diagram_2::Halfedge_around_vertex_circulator vc =
        v->incident_halfedges(),
        vc_begin = vc;
      do {
        Halfedge_handle hh = vc;
        // update the map
        traversed_vertices[Vertex_face_pair(v, hh->face())] = hh;
        // check for reult
        if (hh->face() == f)
          result = hh;
        ++vc;
      } while (vc != vc_begin);
    }
    else
    {
      // take it from the map
      result = iter->second;
    }
    CGAL_assertion(result != Halfedge_handle());
    return result;
  }

  // update halfedge-target flags of new_h that is created on edge on_edge
  // and target-face flags
  // id is the source diagram where on_edge comes from
  // (i.e. the id of the aux information to update)
  void update_halfedge_flags_on_edge(Halfedge_handle new_h, Halfedge_handle on_edge, unsigned int id)
  {
    if(new_h->target()->is_at_open_boundary())
      return;
    Vertex_handle vh;
    Halfedge_handle hh;
    const Object& trg_src = new_h->target()->get_aux_source(id);
    if (assign(vh, trg_src))
    {
          // vh is the target of on_edge, and we can copy the halfedge-target information
      // from on_edge
      copy_halfedge_target_info(on_edge, new_h, id);
            new_h->set_has_equal_aux_data_in_target_and_face
                         (id, on_edge->get_has_equal_data_in_target_and_face());
    }
    else if (assign(hh, trg_src))
    {
      // hh is the "HEMSHECH" of on_edge, so we need to set halfedge_target
      // information to true
      set_halfedge_target_info(new_h, id, true);
            // and target-face information using the original halfedge-face information
            new_h->set_has_equal_aux_data_in_target_and_face
                         (id, on_edge->get_has_equal_data_in_face());
    }
    else
      // this cannot happen, since we need to touch an edge
      CGAL_assertion(false);
  }

  // update halfedge-target flags of new_h that is created inside face in_face
  // and target-face information
  // id is the source diagram where in_face comes from
  // (i.e. the id of the aux information to update)
  void update_halfedge_flags_in_face(Halfedge_handle new_h, Face_handle in_face, unsigned int id)
  {
    if(new_h->target()->is_at_open_boundary())
      return;
    Vertex_handle vh;
    Halfedge_handle hh;
    Face_handle fh;
    // update target
    const Object& trg_src = new_h->target()->get_aux_source(id);
    if (assign(vh, trg_src))
    {
      if (vh->is_isolated())
      {
        copy_halfedge_target_info_from_vertex_face_info(vh, new_h, id);
                    // the target-face information is taken from vertex-face information too
        new_h->set_has_equal_aux_data_in_target_and_face
                                     (id, vh->get_has_equal_data_in_face());
      }
      else
      {
        // we have a vertex vh on the boundary of the face in_face
        // todo: get rid of this calculations: (using unknown value for has_equal flag)
        /*CGAL_assertion_code(
          bool calc_is_equal = vh->is_equal_data(in_face->begin_data(), in_face->end_data());
        )*/
        //bool calc_has_equal = vh->has_equal_data(in_face->begin_data(), in_face->end_data());

        // find the halfedge with target vh on the boundary of in_face
        Halfedge_handle h_of_vh_and_in_face =
               find_halfedge_by_vertex_and_face(vh, in_face);
        // is_equal relationship is easy:
        bool is_equal = h_of_vh_and_in_face->get_is_equal_data_in_face() &&
                        h_of_vh_and_in_face->get_is_equal_data_in_target();
        //CGAL_assertion(is_equal == calc_is_equal);

        // has_equal relationship is problematic in one case:
        bool has_equal =
          h_of_vh_and_in_face->get_has_equal_data_in_target_and_face();

        /*CGAL_assertion(has_equal == calc_has_equal);
        if(has_equal != calc_has_equal)
          return;*/

                    // update halfedge-target flags
        new_h->set_is_equal_aux_data_in_target(id, is_equal);
        new_h->set_has_equal_aux_data_in_target(id, has_equal);
                    // update target-face flag
        new_h->set_has_equal_aux_data_in_target_and_face(id, has_equal);
      }
    }
    else if (assign(hh, trg_src))
      // we should find the halfedge (hh or hh>twin()) that points to face
      // in_face, and check the halfedge-face flags there
    {
      CGAL_assertion(hh->face() == in_face || hh->twin()->face() == in_face);
      if (hh->face() == in_face)
        copy_halfedge_target_info_from_halfedge_face_info(hh, new_h, id);
      else
        copy_halfedge_target_info_from_halfedge_face_info(hh->twin(), new_h, id);
    }
    else
    {
      CGAL_assertion_code(bool b =)
      assign(fh, trg_src);
      CGAL_assertion(b);
      // the edge and target are in the same face, so we set halfedge-target
      // is_equal information to true, and has_equal information according to
      // the face data
      CGAL_assertion(fh == in_face);
      new_h->set_is_equal_aux_data_in_target(id, true);
      new_h->set_has_equal_aux_data_in_target(id, !fh->has_no_data());
      new_h->set_has_equal_aux_data_in_target_and_face(id, !fh->has_no_data());
    }
  }

  Minimization_diagram_2& m_1;
  Minimization_diagram_2& m_2;
  Minimization_diagram_2& m_result;
  Boundary_cache          traversed_vertices;
};

} //namespace CGAL

#endif
