// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ron Wein   <wein@post.tau.ac.il>
//            Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H
#define CGAL_ENVELOPE_DIVIDE_AND_CONQUER_2_H

#include <CGAL/license/Envelope_2.h>

#include <memory>
#include <optional>
#include <vector>
#include <variant>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>

namespace CGAL {

/*! \class Envelope_divide_and_conquer_2
 * A class implementing the divide-and-conquer algorithm for computing the
 * lower (or upper) envelope of a set of curves.
 */
template <typename Traits_, typename Diagram_>
class Envelope_divide_and_conquer_2 {
public:
  using Traits_2 = Traits_;
  using Point_2 = typename Traits_2::Point_2;
  using X_monotone_curve_2 = typename Traits_2::X_monotone_curve_2;
  using Curve_2 = typename Traits_2::Curve_2;

  using Envelope_diagram_1 = Diagram_;

protected:
  using Self = Envelope_divide_and_conquer_2<Traits_2, Envelope_diagram_1>;

  enum Envelope_type {
    LOWER,
    UPPER
  };

  using Vertex_const_handle = typename Envelope_diagram_1::Vertex_const_handle;
  using Vertex_handle = typename Envelope_diagram_1::Vertex_handle;
  using Edge_const_handle = typename Envelope_diagram_1::Edge_const_handle;
  using Edge_handle = typename Envelope_diagram_1::Edge_handle;

  using Curve_pointer_vector = std::vector<X_monotone_curve_2 *>;
  using Curve_pointer_iterator = typename Curve_pointer_vector::iterator;

  using Traits_adaptor_2 = Arr_traits_adaptor_2<Traits_2>;

  // Data members:
  std::shared_ptr<const Traits_2> traits_ptr;
  const Traits_adaptor_2* traits;
  bool own_traits;
  Envelope_type env_type;

  Envelope_divide_and_conquer_2(const Self&) = delete;
  Self& operator=(const Self&) = delete;

public:
  /*! Default constructor. */
  Envelope_divide_and_conquer_2() :
    traits_ptr(std::make_shared<Traits_2>()),
    own_traits(true),
    env_type(LOWER)
  { traits = new Traits_adaptor_2(*traits_ptr); }

  /*! Constructor with a traits object pointer. */
  Envelope_divide_and_conquer_2(const Traits_2* _traits) :
    traits_ptr(_traits ? std::shared_ptr<const Traits_2>(_traits, [](const Traits_2*){}) : nullptr),
    own_traits(false),
    env_type(LOWER)
  { traits = static_cast<const Traits_adaptor_2*>(_traits); }

  /*! Destructor. */
  ~Envelope_divide_and_conquer_2()
  { if (own_traits) delete traits; }

  /*! Construct the lower (or upper) envelope to the given range of curves. */
  template <typename CurvesIterator>
  void insert_curves(CurvesIterator begin, CurvesIterator end, bool type, Envelope_diagram_1& diagram) {
    std::list<std::variant<Point_2, X_monotone_curve_2>> objects;
    std::list<X_monotone_curve_2> x_curves;

    for (auto it = begin; it != end; ++it) {
      objects.clear();
      traits->make_x_monotone_2_object()(*it, std::back_inserter(objects));

      for (auto obj_it = objects.begin(); obj_it != objects.end(); ++obj_it) {
        if (const auto* xcv_ptr = std::get_if<X_monotone_curve_2>(&(*obj_it))) x_curves.push_back(*xcv_ptr);
      }
    }

    insert_x_monotone_curves(x_curves.begin(), x_curves.end(), type, diagram);
  }

  /*! Construct the lower (or upper) envelope to the given range of x-monotone curves. */
  template <typename XCurvesIterator>
  void insert_x_monotone_curves(XCurvesIterator begin, XCurvesIterator end, bool type, Envelope_diagram_1& diagram) {
    env_type = (type ? LOWER : UPPER);

    auto is_vertical = traits->is_vertical_2_object();

    Curve_pointer_vector reg_vec;
    Curve_pointer_vector vert_vec;

    for (auto iter = begin; iter != end; ++iter) {
      if (is_vertical(*iter)) vert_vec.push_back(&(*iter));
      else reg_vec.push_back(&(*iter));
    }

    _construct_envelope_non_vertical(reg_vec.begin(), reg_vec.end(), diagram);

    if (!vert_vec.empty()) _merge_vertical_segments(vert_vec, diagram);
  }

  /*! Get pointer to traits. */
  const Traits_2* get_traits() const { return traits_ptr.get(); }

protected:
  void _construct_envelope_non_vertical(Curve_pointer_iterator begin,
                                        Curve_pointer_iterator end,
                                        Envelope_diagram_1& out_d);

  void _construct_singleton_diagram(const X_monotone_curve_2& cv,
                                    Envelope_diagram_1& out_d);

  void _merge_envelopes(const Envelope_diagram_1& d1,
                        const Envelope_diagram_1& d2,
                        Envelope_diagram_1& out_d);

  Comparison_result _compare_vertices(const Envelope_diagram_1& d1,
                                      const Envelope_diagram_1& d2,
                                      Vertex_const_handle v1,
                                      Vertex_const_handle v2,
                                      bool& same_x) const;

  void _merge_single_interval(Edge_const_handle e,
                              Edge_const_handle other_edge,
                              Vertex_const_handle v, bool v_exists,
                              Comparison_result origin_of_v,
                              const Envelope_diagram_1& in_d,
                              const Envelope_diagram_1& other_d,
                              Envelope_diagram_1& out_d);

  Comparison_result compare_y_at_end(const X_monotone_curve_2& xcv1,
                                     const X_monotone_curve_2& xcv2,
                                     Arr_curve_end curve_end) const;

  void _merge_two_intervals(Edge_const_handle e1, bool is_leftmost1,
                            Edge_const_handle e2, bool is_leftmost2,
                            Vertex_const_handle v, bool v_exists,
                            Comparison_result origin_of_v,
                            const Envelope_diagram_1& d1,
                            const Envelope_diagram_1& d2,
                            Envelope_diagram_1& out_d);

  Vertex_handle _append_vertex(Envelope_diagram_1& diag,
                               const Point_2& p, Edge_const_handle e,
                               const Envelope_diagram_1& in_d);

  class Less_vertical_segment {
  private:
    typename Traits_2::Compare_x_2 comp_x;
    typename Traits_2::Construct_min_vertex_2 min_vertex;

  public:
    Less_vertical_segment(const Traits_2 *traits) :
      comp_x(traits->compare_x_2_object()),
      min_vertex(traits->construct_min_vertex_2_object())
    {}

    bool operator()(const X_monotone_curve_2* cv1, const X_monotone_curve_2* cv2) const
    { return(comp_x(min_vertex(*cv1), min_vertex(*cv2)) == SMALLER); }
  };

  void _merge_vertical_segments(Curve_pointer_vector& vert_vec, Envelope_diagram_1& out_d);

  Vertex_handle _split_edge(Envelope_diagram_1& diag, const Point_2& p, Edge_handle e);
};

} // namespace CGAL

#include <CGAL/Envelope_2/Envelope_divide_and_conquer_2_impl.h>

#endif
