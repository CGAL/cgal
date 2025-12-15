// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Efi Fogel        <efifogel@gmail.com>

#ifndef CGAL_GSP_DO_INTERSECT_AGG_OP_VISITOR_H
#define CGAL_GSP_DO_INTERSECT_AGG_OP_VISITOR_H

#include <vector>

#include <CGAL/license/Boolean_set_operations_2.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op_visitor.h>
#include <CGAL/Default.h>

namespace CGAL {

template <typename Helper_, typename Arrangement_, typename Visitor_ = Default>
class Gps_do_intersect_agg_op_visitor :
  public Gps_agg_op_visitor<
    Helper_, Arrangement_,
    typename Default::Get<Visitor_, Gps_do_intersect_agg_op_visitor<Helper_, Arrangement_, Visitor_>>::type> {
public:
  using Helper = Helper_;
  using Arrangement_2 = Arrangement_;
  using Geometry_traits_2 = typename Helper::Geometry_traits_2;
  using Event = typename Helper::Event;
  using Subcurve = typename Helper::Subcurve;

private:
  using Gt2 = Geometry_traits_2;
  using Arr = Arrangement_2;
  using Self = Gps_do_intersect_agg_op_visitor<Helper, Arr, Visitor_>;
  using Visitor = typename Default::Get<Visitor_, Self>::type;
  using Base = Gps_agg_op_visitor<Helper, Arr, Visitor>;

protected:
  bool m_found_x;

public:
  using Edges_hash = typename Base::Edges_hash;
  using Vertex_handle = typename Base::Vertex_handle;
  using Status_line_iterator = typename Base::Status_line_iterator;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Point_2 = typename Base::Point_2;
  using Multiplicity = typename Base::Multiplicity;

  Gps_do_intersect_agg_op_visitor(Arr* arr, Edges_hash* hash,
                                  std::vector<Vertex_handle>* vertices_vec) :
    Base(arr, hash, vertices_vec),
    m_found_x(false)
  {}

  /*! Update an event that corresponds to a curve endpoint. */
  void update_event(Event* e, const Point_2& end_point, const X_monotone_curve_2& cv, Arr_curve_end cv_end, bool is_new)
  { Base::update_event(e, end_point, cv, cv_end, is_new); }

  /*! Update an event that corresponds to a curve endpoint */
  void update_event(Event* e, const X_monotone_curve_2& cv, Arr_curve_end cv_end, bool is_new )
  { Base::update_event(e, cv, cv_end, is_new); }

  /*! Update an event that corresponds to a curve endpoint */
  void update_event(Event* e, const Point_2& p, bool is_new)
  { Base::update_event(e, p, is_new); }

  /*! Update an event that corresponds to an intersection */
  void update_event(Event* e, Subcurve* sc) { Base::update_event(e, sc); }

  /*! Update an event that corresponds to an intersection between curves */
  void update_event(Event* e, Subcurve* sc1, Subcurve* sc2, bool is_new, Multiplicity multiplicity) {
    if ((multiplicity % 2) == 1) m_found_x = true;
    Base::update_event(e, sc1, sc2, is_new, multiplicity);
  }

  //!
  bool after_handle_event(Event* e, Status_line_iterator iter, bool flag) {
    auto res = Base::after_handle_event(e, iter, flag);
    if (m_found_x) this->surface_sweep()->stop_sweep();
    return res;
  }

  /*! Getter */
  bool found_intersection() { return m_found_x; }
};

} // namespace CGAL

#endif
