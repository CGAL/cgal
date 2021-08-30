// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_FACE_BASE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_FACE_BASE_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Triangulation_ds_face_base_2.h>


namespace CGAL {

template < typename Gt, typename Fb = Triangulation_ds_face_base_2<> >
class Segment_Delaunay_graph_face_base_2
  : public Fb
{
private:
  typedef typename Fb::Triangulation_data_structure   D_S;
  typedef Fb                                          Base;

public:
  // TYPES
  //------
  typedef Gt                                          Geom_traits;
  typedef typename Fb::Vertex_handle                  Vertex_handle;
  typedef typename Fb::Face_handle                    Face_handle;
  typedef D_S                                         Data_structure;
  typedef typename Data_structure::Edge               Edge;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other  Fb2;
    typedef Segment_Delaunay_graph_face_base_2<Gt, Fb2>    Other;
  };

private:
  class Face_data
  {
    bool in_conflict_;
    Sign incircle_sign_;
    // below are the data for the in-place edge list
    bool selected[3];
    Edge next_[3], prev_[3];

    static int sentinel_index() { return -1; }

    static Edge sentinel_edge() {
      return Edge(Face_handle(), sentinel_index());
    }

    // method to initialize the in-place edge list
    void initialize_in_place_edge_list()
    {
      for (int i = 0; i < 3; ++i) {
        selected[i] = false;
        next_[i] = sentinel_edge();
        prev_[i] = sentinel_edge();
      }
    }

  public:
    Face_data() : in_conflict_(false), incircle_sign_(ZERO)
    {
      initialize_in_place_edge_list();
    }

    inline void clear()                    { in_conflict_ = false; }
    inline void mark_in_conflict()         { in_conflict_ = true; }
    inline void set_incircle_sign(Sign s)  { incircle_sign_ = s; }

    inline bool is_clear()       const { return !in_conflict_; }
    inline bool is_in_conflict() const { return in_conflict_; }
    inline Sign incircle_sign()  const { return incircle_sign_; }

    // in-place edge list stuff -- start
    inline Edge next(unsigned int i)     const { return next_[i]; }
    inline Edge previous(unsigned int i) const { return prev_[i]; }

    inline void set_next(unsigned int i, const Edge& next)
    {
      next_[i] = next;
    }

    inline void set_previous(unsigned int i, const Edge& prev)
    {
      prev_[i] = prev;
    }

    inline bool is_in_list(unsigned int i) const {
      return next_[i].second != -1 || prev_[i].second != -1;
    }

    inline void mark_selected(unsigned int i) { selected[i] = true; }
    inline void mark_unselected(unsigned int i) { selected[i] = false; }
    inline bool is_selected(unsigned int i) const { return selected[i]; }
    // in-place edge list stuff -- end
  };


private:
  // MEMBER DATA
  Face_data      _tds_data;

public:
  Segment_Delaunay_graph_face_base_2()
    : Fb() {}

  Segment_Delaunay_graph_face_base_2(Vertex_handle v0,
                                     Vertex_handle v1,
                                     Vertex_handle v2)
    : Fb(v0,v1,v2) {}

  Segment_Delaunay_graph_face_base_2(Vertex_handle v0,
                                     Vertex_handle v1,
                                     Vertex_handle v2,
                                     Face_handle n0,
                                     Face_handle n1,
                                     Face_handle n2)
    : Fb(v0,v1,v2,n0,n1,n2) {}


  typedef Face_data   TDS_data;

  // TDS internal data access functions.
        TDS_data& tds_data()       { return _tds_data; }
  const TDS_data& tds_data() const { return _tds_data; }
};


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_SEGMENT_DELAUNAY_GRAPH_FACE_BASE_2_H
