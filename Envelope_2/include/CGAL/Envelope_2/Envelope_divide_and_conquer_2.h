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

#ifndef CGAL_ENVELOPE_2_ENVELOPE_DIVIDE_AND_CONQUER_2_H
#define CGAL_ENVELOPE_2_ENVELOPE_DIVIDE_AND_CONQUER_2_H

#include <CGAL/license/Envelope_2.h>

#include <memory>
#include <optional>
#include <vector>
#include <variant>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_tags.h>

namespace CGAL {
namespace Envelope_2 {

#define CGAL_VALUE_BASED_POOL 1

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

  enum Envelope_type { LOWER, UPPER };

protected:
  using Self = Envelope_divide_and_conquer_2<Traits_2, Envelope_diagram_1>;

  using Vertex_const_descriptor = typename Envelope_diagram_1::Vertex_const_descriptor;
  using Vertex_descriptor = typename Envelope_diagram_1::Vertex_descriptor;
  using Edge_const_descriptor = typename Envelope_diagram_1::Edge_const_descriptor;
  using Edge_descriptor = typename Envelope_diagram_1::Edge_descriptor;

  using Curve_pointer_vector = std::vector<X_monotone_curve_2*>;
  using Curve_pointer_iterator = typename Curve_pointer_vector::iterator;

  using Traits_adaptor_2 = Arr_traits_adaptor_2<Traits_2>;

  // All sides
  using Left_side_category = typename CGAL::internal::Arr_complete_left_side_category<Traits_2>::Category;
  using Bottom_side_category = typename CGAL::internal::Arr_complete_bottom_side_category<Traits_2>::Category;
  using Top_side_category = typename CGAL::internal::Arr_complete_top_side_category<Traits_2>::Category;
  using Right_side_category = typename CGAL::internal::Arr_complete_right_side_category<Traits_2>::Category;
  using Are_all_sides_oblivious_category = typename Arr_all_sides_oblivious_category<Left_side_category,
                                                                                     Bottom_side_category,
                                                                                     Top_side_category,
                                                                                     Right_side_category>::result;

private:
  // Data members:
  const Traits_adaptor_2* m_traits;     // the traits object.
  bool m_own_traits;                    // whether we own the traits object.
  Envelope_type m_env_type;             // either LOWER or UPPER.

#if CGAL_VALUE_BASED_POOL==1
  std::vector<Envelope_diagram_1> m_diagram_pool;
#else
  std::vector<Envelope_diagram_1*> m_diagram_pool;
#endif

  // copy constructor and assignment operator - not supported.
  Envelope_divide_and_conquer_2(const Self&) = delete;
  Self& operator=(const Self&) = delete;

public:
  /*! default constructor.
   */
  Envelope_divide_and_conquer_2(Envelope_type type = LOWER) :
    m_own_traits(true),
    m_env_type(type)
  { m_traits = new Traits_adaptor_2; }

  /*! constructor with a traits object.
   * \param traits The traits object.
   */
  Envelope_divide_and_conquer_2(const Traits_2* traits) :
    m_own_traits(false),
    m_env_type(LOWER)
  { m_traits = static_cast<const Traits_adaptor_2*>(traits); }

  /*! constructor with a traits object.
   * \param traits The traits object.
   */
  Envelope_divide_and_conquer_2(Envelope_type type, const Traits_2* traits) :
    m_own_traits(false),
    m_env_type(type)
  { m_traits = static_cast<const Traits_adaptor_2*>(traits); }

  /*! destructor.
   */
  ~Envelope_divide_and_conquer_2() { if (m_own_traits) delete m_traits; }

  /*! constructs the lower (or upper) envelope to the given range of curves.
   * \param begin An iterator pointing at the beginning of the curves range.
   * \param end A past-the-end iterator for the curves range.
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <typename CurvesIterator>
  void insert_curves(CurvesIterator begin, CurvesIterator end, Envelope_diagram_1& diagram) {
    // Subdivide the curves into x-monotone subcurves.
    std::list<std::variant<Point_2, X_monotone_curve_2>> objects;
    std::list<X_monotone_curve_2> x_curves;

    for (auto it = begin; it != end; ++it) {
      // Split the current curve to x-monotone subcurves.
      objects.clear();
      m_traits->make_x_monotone_2_object()(*it, std::back_inserter(objects));

      for (auto obj_it = objects.begin(); obj_it != objects.end(); ++obj_it) {
        if (const auto* xcv_ptr = std::get_if<X_monotone_curve_2>(&(*obj_it))) x_curves.push_back(*xcv_ptr);
      }
    }

    // Construct the envelope of the x-monotone curves.
    insert_x_monotone_curves(x_curves.begin(), x_curves.end(), diagram);
  }

  /*! constructs the lower (or upper) envelope to the given range of
   * x-monotone curves.
   * \param begin An iterator pointing at the beginning of the curves range.
   * \param end A past-the-end iterator for the curves range.
   * \param diagram Output: The minimization (or maximization) diagram.
   */
  template <typename XCurvesIterator>
  void insert_x_monotone_curves(XCurvesIterator begin, XCurvesIterator end, Envelope_diagram_1& diagram) {
    // Separate the regular curves from the vertical ones.
    auto is_vertical = m_traits->is_vertical_2_object();

    Curve_pointer_vector reg_vec;
    Curve_pointer_vector vert_vec;

    // Pre-reserve to avoid vector buffer growth
    auto num_curves = std::distance(begin, end);
    reg_vec.reserve(num_curves);
    vert_vec.reserve(num_curves);

    for (auto iter = begin; iter != end; ++iter) {
      if (is_vertical(*iter)) vert_vec.push_back(&(*iter));
      else reg_vec.push_back(&(*iter));
    }

    // Construct the envelope for the non-vertical curves.
    _construct_envelope_non_vertical(reg_vec.begin(), reg_vec.end(), diagram);

    // Merge the vertical segments.
    if (! vert_vec.empty()) _merge_vertical_segments(vert_vec, diagram);
  }

  /*! obtains the traits object.
   * \return A pointer to the traits object.
   */
  const Traits_2* geometry_traits() const { return m_traits; }

protected:
  /*! constructs the lower/upper envelope of the given list of non-vertical curves.
   * \param begin The first x-monotone curve.
   * \param end A past-the-end iterator for the curves.
   * \param out_d Output: The minimization (or maximization) diagram.
   */
  void _construct_envelope_non_vertical(Curve_pointer_iterator begin, Curve_pointer_iterator end,
                                        Envelope_diagram_1& out_d);

  /*!
   */
#if CGAL_VALUE_BASED_POOL==1
  void _construct_envelope_non_vertical_pooled(Curve_pointer_iterator begin, Curve_pointer_iterator end,
                                               Envelope_diagram_1& out_d, std::size_t& pool_active_count);
#else
  void _construct_envelope_non_vertical_pooled(Curve_pointer_iterator begin, Curve_pointer_iterator end,
                                               Envelope_diagram_1& out_d);
#endif

  /*! constructs a diagram of a single curve.
   * \param cv The x-monotone curve.
   * \param out_d Output: The minimization (or maximization) diagram.
   * This is the implementation for the case where all 4 boundary sides are oblivious.
   */
  void _construct_singleton_diagram(const X_monotone_curve_2& xcv, Envelope_diagram_1& out_d,
                                    Arr_all_sides_oblivious_tag);

  /*! constructs a diagram of a single curve.
   * \param cv The x-monotone curve.
   * \param out_d Output: The minimization (or maximization) diagram.
   * This is the implementation for the case where at least one of the 4 boundary sides are not oblivious.
   */
  void _construct_singleton_diagram(const X_monotone_curve_2& xcv, Envelope_diagram_1& out_d,
                                    Arr_not_all_sides_oblivious_tag);

  /* merges two minimization (or maximization) diagrams.
   * \param d1 The first diagram,
   *           representing the envelope of the curve set C1.
   * \param d2 The second diagram,
   *           representing the envelope of the curve set C2.
   * \param out_d Output: The merged diagram, representing the envelope of
   *                      the union of C1 and C2.
   */
  void _merge_envelopes(const Envelope_diagram_1& d1, const Envelope_diagram_1& d2, Envelope_diagram_1& out_d);

  /*! compares two vertices.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \param same_x Output parameter: TRUE iff x(v1) = x(v2).
   * \return SMALLER if x(v1) < x(v2). Or, in case x(v1) = x(v2), and
   *                 - we compute the lower envelope, and y(v1) < y(v2),
   *                 - we compute the upper envelope, and y(v1) > y(v2).
   *         LARGER if x(v1) > x(v2). Or, in case x(v1) = x(v2), and
   *                - we compute the lower envelope, and y(v1) > y(v2),
   *                - we compute the upper envelope, and y(v1) < y(v2).
   *         EQUAL if v1 = v2.
   */
  Comparison_result _compare_vertices(const Envelope_diagram_1& d1, const Envelope_diagram_1& d2,
                                      Vertex_const_descriptor v1, Vertex_const_descriptor v2, bool& same_x) const;

  /*! deals with an interval which is non-empty in one of the merged diagrams
   * and empty in the other.
   * \param e The non-empty edge.
   * \param other_edge The empty edge.
   * \param v The next vertex to the right.
   * \param v_exists Whether the next vertex exists.
   * \param origin_of_v The origin of v: SMALLER if it is from e,
   *                    LARGER if it is from other_edge.
   *                    EQUAL result means that both edges have vertex at
   *                    the same place.
   * \param out_d The merged diagram.
   */
  void _merge_single_interval(Edge_const_descriptor e, Edge_const_descriptor other_edge,
                              Vertex_const_descriptor v, bool v_exists, Comparison_result origin_of_v,
                              const Envelope_diagram_1& in_d, const Envelope_diagram_1& other_d,
                              Envelope_diagram_1& out_d);

  /*! compares the \f$y\f$-coordinates of two curves at their endpoints
   * The function compares the \f$y\f$ values of two curves with a joint
   * range of \f$x\f$ values, at the end of the joint range.
   * \param xcv1 The first curve
   * \param xcv2 The second curve
   * \param curve_end `ARR_MIN_END` - compare the \f$y\f$ value of the smaller endpoint,
   *                  `ARR_MAX_END` - compare the \f$y\f$ value of the larger endpoint.
   * \pre The two \f$x\f$-monotone curves need to have a partially overlapping \f$x\f$-ranges.
   * \return `SMALLER` - the end of `xcv1` is below the end of `xcv2`.
   *         `EQUAL`   - the end of `xcv1` and the end of `xcv2` have the same \f$y\f$ coordinates.
   * \       `LARGER`  - the end of `xcv1` is above the end of `xcv2`.
   */
  Comparison_result compare_y_at_end(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                     Arr_curve_end curve_end, Arr_all_sides_oblivious_tag) const;
  Comparison_result compare_y_at_end(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                     Arr_curve_end curve_end, Arr_not_all_sides_oblivious_tag) const;

  /*! merges two non-empty intervals into the merged diagram.
   * \param e1 The first non-empty edge.
   * \param is_leftmost1 Is it the leftmost edge in its diagram.
   * \param e2 The second non-empty edge.
   * \param is_leftmost2 Is it the leftmost edge in its diagram.
   * \param v The next vertex.
   * \param v_exists Whether such a vertex exists.
   * \param origin_of_v The origin of v: SMALLER if it is from e1,
   *                    otherwise it is from e2. EQUAL result means that
   *                    both diagram have vertex at the same place (but v
   *                    is still taken from e2.
   * \param out_d The merged diagram.
   */
  void _merge_two_intervals(Edge_const_descriptor e1, bool is_leftmost1, Edge_const_descriptor e2, bool is_leftmost2,
                            Vertex_const_descriptor v, bool v_exists, Comparison_result origin_of_v,
                            const Envelope_diagram_1& d1, const Envelope_diagram_1& d2, Envelope_diagram_1& out_d);

  /*! appends a vertex to the given diagram: The new vertex that represents the
   * given point as the new rightmost vertex of the diagram. The edge
   * between the current rightmost vertex and the new one contains the same
   * curves as the input edge.
   * \param diag The diagram.
   * \param p The point that the new vertex is associated with.
   * \param e The input edge.
   * \return A descriptor for the vertex.
   */
  Vertex_descriptor _append_vertex(Envelope_diagram_1& diag, const Point_2& p, Edge_const_descriptor e,
                                   const Envelope_diagram_1& in_d);

  /*! \struct
   * A functor used to sort vertical segments by their x-coordinate.
   */
  class Less_vertical_segment {
  private:
    typename Traits_2::Compare_x_2 m_comp_x;
    typename Traits_2::Construct_min_vertex_2 m_min_vertex;

  public:
    Less_vertical_segment(const Traits_2* traits) :
      m_comp_x(traits->compare_x_2_object()),
      m_min_vertex(traits->construct_min_vertex_2_object())
    {}

    bool operator()(const X_monotone_curve_2* cv1, const X_monotone_curve_2* cv2) const
    { return(m_comp_x(m_min_vertex(*cv1), m_min_vertex(*cv2)) == SMALLER); }
  };

  /*! merges the vertical segments into the lower/upper envelope given as a
   * minimization (or maximization) diagram.
   * \param vert_vec The list of vertical segments.
   * \param out_d The input minimization (or maximization) diagram.
   *              The function merges the vertical segments into this diagram.
   */
  void _merge_vertical_segments(Curve_pointer_vector& vert_vec, Envelope_diagram_1& out_d);

  /*! splits a given diagram edge by inserting a vertex in its interior.
   * \param diag The diagram.
   * \param p The point that the new vertex is associated with.
   * \param e The edge to split.
   * \return A descriptor for the vertex.
   */
  Vertex_descriptor _split_edge(Envelope_diagram_1& diag, const Point_2& p, Edge_descriptor e);
};

} // namespace Envelope_2
} // namespace CGAL

#include <CGAL/Envelope_2/Envelope_divide_and_conquer_2_impl.h>

#endif
