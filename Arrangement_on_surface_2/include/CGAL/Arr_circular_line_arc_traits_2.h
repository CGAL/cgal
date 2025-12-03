// Copyright (c) 2003,2004,2005,2006,2007,2008,2009,2010,2011 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473
// (ECG - Effective Computational Geometry for Curves and Surfaces)
// and a STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_VARIANT_TRAITS_2_H
#define CGAL_CIRCULAR_KERNEL_VARIANT_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * This file was developed at Inria, France, and copied over to the
 * Arrangement_2 package, which it is now part of. It contains a traits
 * class for the arrangement package that handles circular and linear curves.
 * It is based on the circular kernel.
 *
 * \todo Fix the circular-kernel make-x-monotone functor to use modern variant
 *       instead of the legacy CGAL::Object. Then, eliminate the special
 *       implementation here and directly use the kernel functor instead.
 */

#include <vector>

#include <variant>

#include <CGAL/basic.h>
#include <CGAL/Arr_tags.h>

namespace CGAL {

namespace VariantFunctors {

// Takes an iterator range of Object(Line/Circular_arc/Point),
// returns a variant of Line, Circular_arc, and Point_2.
template <typename CK, typename Arc1, typename Arc2, typename OutputIterator>
OutputIterator object_to_object_variant(const std::vector<CGAL::Object>& res1,
                                        OutputIterator res2) {
  using Point_2 = typename CK::Circular_arc_point_2;
  using X_monotone_curve_2 = std::variant<Arc1, Arc2>;
  using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

  for (auto it = res1.begin(); it != res1.end(); ++it) {
    if (const Arc1* arc = CGAL::object_cast<Arc1>(&*it)) {
      std::variant<Arc1, Arc2> v = *arc;
      *res2++ = Make_x_monotone_result(v);
    }
    else if (const Arc2* line = CGAL::object_cast<Arc2>(&*it)) {
      std::variant<Arc1, Arc2> v = *line;
      *res2++ = Make_x_monotone_result(v);
    }
    else if (const Point_2* p = CGAL::object_cast<Point_2>(&*it)) {
      *res2++ = Make_x_monotone_result(*p);
    }
    else CGAL_error();
  }
  return res2;
}

template <typename CircularKernel, typename Arc1, typename Arc2>
class Compare_y_to_right_2 {
public:
  using result_type = CGAL::Comparison_result;
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;

  result_type operator()(const std::variant<Arc1, Arc2>& a1,
                         const std::variant<Arc1, Arc2>& a2,
                         const Circular_arc_point_2& p) const {
    if (const Arc1* arc1 = std::get_if<Arc1>(&a1)) {
      if (const Arc1* arc2 = std::get_if<Arc1>(&a2)) {
        return CircularKernel().compare_y_to_right_2_object()(*arc1, *arc2, p);
      }
      else {
        const Arc2* arc2e = std::get_if<Arc2>(&a2);
        return CircularKernel().compare_y_to_right_2_object()(*arc1, *arc2e, p);
      }
    }
    const Arc2* arc1 = std::get_if<Arc2>(&a1);
    if (const Arc1* arc2 = std::get_if<Arc1>(&a2)) {
      return CircularKernel().compare_y_to_right_2_object()(*arc1, *arc2, p);
    }
    const Arc2* arc2e = std::get_if<Arc2>(&a2);
    return CircularKernel().compare_y_to_right_2_object()(*arc1, *arc2e, p);
  }
};

template <typename CircularKernel>
class Variant_Equal_2 {
public:
  template <typename T>
  bool operator()(const T& a0, const T& a1) const
  { return CircularKernel().equal_2_object()(a0,a1); }

  template <typename T1, typename T2>
  bool operator()(const T1& , const T2&) const
  { return false; }
};

template <typename CircularKernel, class Arc1, class Arc2>
class Equal_2 : public CircularKernel::Equal_2 {
public:
  using Curve_2 = std::variant< Arc1, Arc2>;
  using result_type = bool;
  using CircularKernel::Equal_2::operator();
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
  using Line_arc_2 = typename CircularKernel::Line_arc_2;
  using Circular_arc_2 = typename CircularKernel::Circular_arc_2;
  using CK_Equal_2 = typename CircularKernel::Equal_2;

  result_type operator()(const Circular_arc_point_2& p0,
                         const Circular_arc_point_2& p1) const
  { return CK_Equal_2()(p0, p1); }

  result_type operator()(const Circular_arc_2& a0, const Circular_arc_2& a1) const
  { return CK_Equal_2()(a0, a1); }

  result_type operator()(const Line_arc_2& a0, const Line_arc_2& a1) const
  { return CK_Equal_2()(a0, a1); }

  result_type operator()(const Line_arc_2& /*a0*/, const Circular_arc_2& /*a1*/) const
  { return false; }

  result_type operator()(const Circular_arc_2& /*a0*/, const Line_arc_2& /*a1*/) const
  { return false; }

  result_type operator()(const Curve_2& a0, const Curve_2& a1) const
  { return std::visit(Variant_Equal_2<CircularKernel>(), a0, a1); }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Compare_y_at_x_2 {
public:
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
  using result_type = CGAL::Comparison_result;

  result_type operator()(const Circular_arc_point_2& p,
                         const std::variant< Arc1, Arc2>& A1) const {
    if (const Arc1* arc1 = std::get_if<Arc1>(&A1)){
      return CircularKernel().compare_y_at_x_2_object()(p, *arc1);
    }
    else {
      const Arc2* arc2 = std::get_if<Arc2>(&A1);
      return CircularKernel().compare_y_at_x_2_object()(p, *arc2);
    }
  }
};

template <typename CircularKernel>
class Variant_Do_overlap_2 {
public:
  template <typename T>
  bool operator()(const T& a0, const T& a1) const
  { return CircularKernel().do_overlap_2_object()(a0, a1); }

  template <typename T1, typename T2>
  bool operator()(const T1&, const T2&) const
  { return false; }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Do_overlap_2 {
public:
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
  using result_type = bool;

  result_type operator()(const std::variant< Arc1, Arc2>& A0,
                         const std::variant< Arc1, Arc2>& A1) const
  { return std::visit(Variant_Do_overlap_2<CircularKernel>(), A0, A1); }
};

//! A functor for subdividing curves into x-monotone curves.
template <typename CircularKernel, typename Arc1, typename Arc2>
class Make_x_monotone_2 {
public:
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;

  template <typename OutputIterator, typename Not_X_Monotone>
  OutputIterator operator()(const std::variant<Arc1, Arc2, Not_X_Monotone>& A,
                            OutputIterator res) const {
    if ( const Arc1* arc1 = std::get_if<Arc1>(&A)) {
      return CircularKernel().make_x_monotone_2_object()(*arc1, res);
    }
    else {
      const Arc2* arc2 = std::get_if<Arc2>(&A);
      return CircularKernel().make_x_monotone_2_object()(*arc2, res);
    }
  }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Intersect_2 {
public:
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;

  template <typename OutputIterator>
  OutputIterator operator()(const std::variant< Arc1, Arc2>& c1,
                            const std::variant< Arc1, Arc2>& c2,
                            OutputIterator oi) const {
    if (const Arc1* arc1 = std::get_if<Arc1>(&c1)) {
      if ( const Arc1* arc2 = std::get_if<Arc1>(&c2)) {
        return CircularKernel().intersect_2_object()(*arc1, *arc2, oi);
      }
      const Arc2* arc2 = std::get_if<Arc2>(&c2);
      return CircularKernel().intersect_2_object()(*arc1, *arc2, oi);
    }

    const Arc2* arc1e = std::get_if<Arc2>(&c1);
    if ( const Arc1* arc2 = std::get_if<Arc1>(&c2)) {
      return CircularKernel().intersect_2_object()(*arc1e, *arc2, oi);
    }
    const Arc2* arc2 = std::get_if<Arc2>(&c2);
    return CircularKernel().intersect_2_object()(*arc1e, *arc2, oi);
  }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Split_2 {
public:
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
  using result_type = void;

  result_type operator()(const std::variant< Arc1, Arc2>& A,
                         const Circular_arc_point_2& p,
                         std::variant< Arc1, Arc2>& ca1,
                         std::variant< Arc1, Arc2>& ca2) const {
    // TODO : optimize by extracting the references from the variants ?
    if (const Arc1* arc1 = std::get_if<Arc1>(&A)) {
      Arc1 carc1;
      Arc1 carc2;
      CircularKernel().split_2_object()(*arc1, p, carc1, carc2);
      ca1 = carc1;
      ca2 = carc2;
      return ;
    }
    else {
      const Arc2* arc2 = std::get_if<Arc2>(&A);
      Arc2 cline1;
      Arc2 cline2;
      CircularKernel().split_2_object()(*arc2, p, cline1, cline2);
      ca1 = cline1;
      ca2 = cline2;
      return ;
    }
  }
};

template <typename CircularKernel>
class Variant_Construct_min_vertex_2 {
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
public:
  using result_type = Circular_arc_point_2;
  // using qualified_result_type = const result_type& ;

  template <typename T>
  //std::remove_reference_t<qualified_result_type>
  Circular_arc_point_2 operator()(const T& a) const {
    //CGAL_kernel_precondition(CircularKernel().compare_xy_2_object()(a.left(), a.right())==CGAL::SMALLER);
    return CircularKernel().construct_circular_min_vertex_2_object()(a);
  }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Construct_min_vertex_2 {
  //: public Has_qrt
  using Point_2 = typename CircularKernel::Circular_arc_point_2;

public:
  using result_type = Point_2;
  // using qualified_result_type = const result_type&;

  //std::remove_reference_t<qualified_result_type>
  result_type operator() (const std::variant< Arc1, Arc2>& cv) const
  { return std::visit(Variant_Construct_min_vertex_2<CircularKernel>(), cv); }
};

template <typename CircularKernel>
class Variant_Construct_max_vertex_2 {
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;
public:
  using result_type = Circular_arc_point_2;
  // using qualified_result_type = const result_type&;

  template <typename T>
    //std::remove_reference_t<qualified_result_type>
  Circular_arc_point_2 operator()(const T& a) const {
    //CGAL_kernel_precondition(CircularKernel().compare_xy_2_object()(a.left(), a.right())==CGAL::SMALLER);
    return (CircularKernel().construct_circular_max_vertex_2_object()(a));
  }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Construct_max_vertex_2 {
  //: public Has_qrt
  using Point_2 = typename CircularKernel::Circular_arc_point_2;

public:
  /*! obtains the right endpoint of the x-monotone curve (segment).
   * \param cv The curve.
   * \return The right endpoint.
   */
  using result_type = Point_2;
  // using qualified_result_type = const result_type&;

  //std::remove_reference<qualified_result_type>
  result_type operator() (const std::variant<Arc1, Arc2>& cv) const
  { return std::visit(Variant_Construct_max_vertex_2<CircularKernel>(), cv); }
};

template <typename CircularKernel>
class Variant_Is_vertical_2 {
public:
  template <typename T>
  bool operator()(const T& a) const
  { return CircularKernel().is_vertical_2_object()(a); }
};

template <typename CircularKernel, typename Arc1, typename Arc2>
class Is_vertical_2 {
public:
  using result_type = bool;

  bool operator() (const std::variant<Arc1, Arc2>& cv) const
  { return std::visit(Variant_Is_vertical_2<CircularKernel>(), cv); }
};

}

// an empty class used to have different types between Curve_2 and X_monotone_curve_2
// in Arr_circular_line_arc_traits_2.
namespace internal_Argt_traits {

struct Not_X_Monotone{};

inline std::ostream& operator << (std::ostream& os, const Not_X_Monotone&) { return os; }

}

/// Traits class for CGAL::Arrangement_2 (and similar) based on a CircularKernel.

template <typename CircularKernel>
class Arr_circular_line_arc_traits_2 {
  using Self = Arr_circular_line_arc_traits_2<CircularKernel>;

  using Arc1 = typename CircularKernel::Line_arc_2;
  using Arc2 = typename CircularKernel::Circular_arc_2;

public:
  using Kernel = CircularKernel;
  using Circular_arc_point_2 = typename CircularKernel::Circular_arc_point_2;

  using Point = typename CircularKernel::Circular_arc_point_2;
  using Point_2 = typename CircularKernel::Circular_arc_point_2;

  using Multiplicity = std::size_t;

  using Has_left_category = CGAL::Tag_false;
  using Has_merge_category = CGAL::Tag_false;
  using Has_do_intersect_category = CGAL::Tag_false;

  using Left_side_category = Arr_oblivious_side_tag;
  using Bottom_side_category = Arr_oblivious_side_tag;
  using Top_side_category = Arr_oblivious_side_tag;
  using Right_side_category = Arr_oblivious_side_tag;

  using Not_X_Monotone = internal_Argt_traits::Not_X_Monotone;

  using Curve_2 = std::variant<Arc1, Arc2, Not_X_Monotone>;
  using X_monotone_curve_2 = std::variant<Arc1, Arc2>;

private:
  CircularKernel ck;

public:
  Arr_circular_line_arc_traits_2(const CircularKernel& k = CircularKernel()) : ck(k) {}

  using Compare_x_2 = typename CircularKernel::Compare_x_2;
  using Compare_xy_2 = typename CircularKernel::Compare_xy_2;
  using Construct_min_vertex_2 = typename  VariantFunctors::Construct_min_vertex_2<CircularKernel, Arc1, Arc2>;
  using Construct_max_vertex_2 = VariantFunctors::Construct_max_vertex_2<CircularKernel, Arc1, Arc2>;
  using Is_vertical_2 = VariantFunctors::Is_vertical_2<CircularKernel, Arc1, Arc2>;
  using Compare_y_at_x_2 = VariantFunctors::Compare_y_at_x_2<CircularKernel, Arc1, Arc2>;
  using Compare_y_at_x_right_2 = VariantFunctors::Compare_y_to_right_2<CircularKernel, Arc1, Arc2>;
  using Equal_2 = VariantFunctors::Equal_2<CircularKernel, Arc1, Arc2>;
  using Make_x_monotone_2 = VariantFunctors::Make_x_monotone_2<CircularKernel, Arc1, Arc2>;
  using Split_2 = VariantFunctors::Split_2<CircularKernel, Arc1, Arc2>;
  using Intersect_2 = VariantFunctors::Intersect_2<CircularKernel, Arc1, Arc2>;

  Compare_x_2 compare_x_2_object() const
  { return ck.compare_x_2_object(); }

  Compare_xy_2 compare_xy_2_object() const
  { return ck.compare_xy_2_object(); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(); }

  Equal_2 equal_2_object() const
  { return Equal_2(); }

  Make_x_monotone_2 make_x_monotone_2_object() const
  { return Make_x_monotone_2(); }

  Split_2 split_2_object() const
  { return Split_2(); }

  Intersect_2 intersect_2_object() const
  { return Intersect_2(); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2();}
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
