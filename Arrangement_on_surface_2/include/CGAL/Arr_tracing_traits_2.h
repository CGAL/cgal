// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Shepard Liu    <shepard0liu@gmail.com>

#ifndef CGAL_ARR_TRACING_TRAITS_H
#define CGAL_ARR_TRACING_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 */

#include <iostream>
#include <list>
#include <type_traits>
#include <variant>
#include <memory>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include "CGAL/Arr_has.h"

namespace CGAL {
namespace aos2 {
namespace internal {

/// `Parameter_space_in_x_2`
//@{

//! Fallback selected `Parameter_space_in_x_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_parameter_space_in_x_2 {};

//! Partial specialization selected if `BaseTraits::`Parameter_space_in_x_2 is defined.
template <typename BaseTraits, typename Derived>
class Tracing_parameter_space_in_x_2<BaseTraits, Derived,
                                     std::enable_if_t<has_parameter_space_in_x_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone curve
  * lies on a boundary of the parameter space along the \f$x\f$-axis.
  */
  class Parameter_space_in_x_2 {
  private:
    typename Base::Parameter_space_in_x_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Parameter_space_in_x_2(const Base& base, bool enabled = true) :
      m_object(base.parameter_space_in_x_2_object()), m_enabled(enabled) {}

    /*! operates
    * \param xcv the curve the end of which is tested.
    * \param ce the curve-end identifier.
    * \return the boundary type.
    */
    Arr_parameter_space operator()(const typename Base::X_monotone_curve_2& xcv, Arr_curve_end ce) const {
      if (! m_enabled) return m_object(xcv, ce);
      std::cout << "parameter_space_in_x" << std::endl
                << "  xcv: " << xcv << ", ce: " << ce << std::endl;
      Arr_parameter_space bt = m_object(xcv, ce);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }

    /*! A functor that obtains the parameter space at a point along the
    * \f$x\f$-axis. Every non-interior point is assumed to lie on the
    * left-right identification. Points at the poles additionally lie on the
    * bottom or top boundary.
    * \param p the point.
    * \return the parameter space at `p`.
    */
    Arr_parameter_space operator()(const typename Base::Point_2& p) const {
      if (! m_enabled) return m_object(p);
      std::cout << "parameter_space_in_x" << std::endl
                << "  p: " << p << std::endl;
      Arr_parameter_space bt = m_object(p);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! obtains a `Parameter_space_in_x_2` function object.
   */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Parameter_space_in_x_2(derived->traits(), derived->parameter_space_in_x_op());
  }
};

//@}

/// `Parameter_space_in_y_2`
//@{

//! Fallback selected `Parameter_space_in_y_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_parameter_space_in_y_2 {};

//! Partial specialization selected if `BaseTraits::Parameter_space_in_y_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_parameter_space_in_y_2<BaseTraits, Derived,
                                     std::enable_if_t<has_parameter_space_in_y_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone arc
   * lies on a boundary of the parameter space along the \f$y\f$-axis.
   */
  class Parameter_space_in_y_2 {
  private:
    typename Base::Parameter_space_in_y_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Parameter_space_in_y_2(const Base& base, bool enabled = true) :
      m_object(base.parameter_space_in_y_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv the curve the end of which is tested.
     * \param ce the curve-end identifier.
     * \return the boundary type.
     */
    Arr_parameter_space operator()(const typename Base::X_monotone_curve_2& xcv, Arr_curve_end ce) const {
      if (! m_enabled) return m_object(xcv, ce);
      std::cout << "parameter_space_in_y" << std::endl
                << "  ce: " << ce << ", xcv: " << xcv << std::endl;
      Arr_parameter_space bt = m_object(xcv, ce);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }

    /*! operates
     * \param p the point.
     * \return the boundary type.
     */
    Arr_parameter_space operator()(const typename Base::Point_2& p) const {
      if (! m_enabled) return m_object(p);
      std::cout << "parameter_space_in_y" << std::endl
                << "  point: " << p << std::endl;
      Arr_parameter_space bt = m_object(p);
      std::cout << "  result: " << bt << std::endl;
      return bt;
    }
  };

  /*! obtains a `Parameter_space_in_y_2` function object.
   */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Parameter_space_in_y_2(derived->traits(), derived->parameter_space_in_y_op());
  }
};

//@}

/// `Make_x_monotone_2`
//@{

//! Fallback selected `Make_x_monotone_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_make_x_monotone_2 {};

//! Fallback selected `Make_x_monotone_2` is not defined in the base traits.
template <typename BaseTraits, typename Derived>
class Tracing_make_x_monotone_2<BaseTraits, Derived, std::enable_if_t<has_make_x_monotone_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that subdivides a curve into \f$x\f$-monotone curves.
  class Make_x_monotone_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
    using Curve_2 = typename Base::Curve_2;

  private:
    typename Base::Make_x_monotone_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Make_x_monotone_2(const Base& base, bool enabled = true) :
      m_object(base.make_x_monotone_2_object()), m_enabled(enabled) {}

    /*! subdivides a given curve into \f$x\f$-monotone subcurves and insert them into a given output iterator.
     * \param cv the curve.
     * \param oi an output iterator for the result. Its value type is a variant
     *           that wraps `Point_2` or `X_monotone_curve_2` objects.
     * \return the output iterator.
     */
    template<typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
      if (! m_enabled) return m_object(cv, oi);
      std::cout << "make_x_monotone" << std::endl
                << "  cv: " << cv << std::endl;

      using Make_x_monotone_result = std::variant<Point_2, X_monotone_curve_2>;

      std::list<Make_x_monotone_result> container;
      m_object(cv, std::back_inserter(container));
      if (container.empty()) return oi;

      std::size_t i = 0;
      for (auto it = container.begin(); it != container.end(); ++it) {
        if (const auto* xcv = std::get_if<X_monotone_curve_2>(&*it)) {
          std::cout << "  result[" << i++ << "]: xcv: " << *xcv << std::endl;
          continue;
        }

        if (const auto* p = std::get_if<Point_2>(&*it)) {
          std::cout << "  result[" << i++ << "]: p: " << *p << std::endl;
          continue;
        }

        CGAL_error();
      }

      for (auto it = container.begin(); it != container.end(); ++it) *oi++ = *it;
      container.clear();
      return oi;
    }
  };

  /*! obtains a `Make_x_monotone_2` function object.
   */
  Make_x_monotone_2 make_x_monotone_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Make_x_monotone_2(derived->traits(), derived->make_x_monotone_op());
  }
};

//@}

/// `Split_2`
//@{

//! Fallback selected `Split_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_split_2 {};

//! Partial specialization selected if `BaseTraits::Split_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_split_2<BaseTraits, Derived, std::enable_if_t<has_split_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that splits an \f$x\f$-monotone curve at a point.
  class Split_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Split_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Split_2(const Base& base, bool enabled = true) :
      m_object(base.split_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv the curve to split.
     * \param p the split point.
     * \param xcv1 the left resulting subcurve (`p` is its right endpoint)..
     * \param xcv2 the right resulting subcurve (`p` is its left endpoint)..
     * \pre `p` lies on `cv` but is not one of its end-points.
     */
    void operator()(const X_monotone_curve_2& xcv, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const {
      if (! m_enabled) {
        m_object(xcv, p, xcv1, xcv2);
        return;
      }
      std::cout << "split: " << std::endl
                << "  xcv: " << xcv << std::endl
                << "  p: " << p << std::endl;
      m_object(xcv, p, xcv1, xcv2);
      std::cout << "  result xcv1: " << xcv1 << std::endl
                << "         xcv2: " << xcv2 << std::endl;
    }
  };

  /*! obtains a `Split_2` function object.
   */
  Split_2 split_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Split_2(derived->traits(), derived->split_op());
  }
};

//@}

/// `Do_intersect_2`
//@{

//! Fallback selected `Do_intersect_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_do_intersect_2 {};

//! Partial specialization selected if `BaseTraits::Do_intersect_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_do_intersect_2<BaseTraits, Derived, std::enable_if_t<has_do_intersect_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether two \f$x\f$-monotone curves intersect.
  class Do_intersect_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Do_intersect_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Do_intersect_2(const Base& base, bool enabled = true) :
      m_object(base.do_intersect_2_object()), m_enabled(enabled) {}

    /*! determines whether two given curves intersect.
     * \param xcv1 the first curve.
     * \param xcv2 the ssecond curve.
     * \param consider_common_endpoints indicates whether common endpoints should be counted as intersections.
     * \return `true` if `consider_common_endpoints` is true and `xcv1` and `xcv2` intersect or if
     *  `consider_common_endpoints` is `false and at least one of the interiors of `xcv1` and `xcv2` intersect,
     *   and `false` otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                    bool consider_common_endpoints = true) const {
      if (! m_enabled) return m_object(xcv1, xcv2, consider_common_endpoints);

      std::cout << "do_intersect" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl
                << "  consider_common_endpoints: " << consider_common_endpoints << std::endl;
      auto res = m_object(xcv1, xcv2, consider_common_endpoints);
      std::cout << "  result: " << res << std::endl;
      return res;
    }
  };

  /*! obtains a `Do_intersect_2` function object.
   */
  Do_intersect_2 do_intersect_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Do_intersect_2(derived->traits(), derived->do_intersect_op());
  }
};

//@}

/// `Intersect_2`
//@{

//! Fallback selected `Intersect_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_intersect_2 {};

//! Partial specialization selected if `BaseTraits::Intersect_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_intersect_2<BaseTraits, Derived, std::enable_if_t<has_intersect_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  using Multiplicity = typename Base::Multiplicity;

  //! A functor that computes intersections between two \f$x\f$-monotone curves.
  class Intersect_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Intersect_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Intersect_2(const Base& base, bool enabled = true) : m_object(base.intersect_2_object()), m_enabled(enabled) {}

    /*! computes the intersections of the two given curves and insert them into
     * a given output iterator.
     * \param xcv1 the first curve.
     * \param xcv2 the ssecond curve.
     * \param oi the output iterator for the result. It value type is a variant
     *           that wraps an \f$x\f$-monotone overlapping curve or a pair that
     *           consists of the intersection point and its multiplicity.
     * \return the past-the-end output iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2, OutputIterator oi) const {
      using Intersection_point = std::pair<Point_2, Multiplicity>;
      using Intersection_result =
        std::variant<Intersection_point, X_monotone_curve_2>;

      if (! m_enabled) return m_object(xcv1, xcv2, oi);

      std::cout << "intersect" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      std::list<Intersection_result> container;
      m_object(xcv1, xcv2, std::back_inserter(container));
      if (container.empty()) return oi;

      std::size_t i = 0;
      for (const auto& item : container) {
        const X_monotone_curve_2* xcv = std::get_if<X_monotone_curve_2>(&item);
        if (xcv != nullptr) {
          std::cout << "  result[" << i++ << "]: xcv: " << *xcv << std::endl;
          *oi++ = *xcv;
          continue;
        }

        const Intersection_point* ip = std::get_if<Intersection_point>(&item);
        if (ip != nullptr) {
          std::cout << "  result[" << i++ << "]: p: " << ip->first << ", multiplicity: " << ip->second << std::endl;
          *oi++ = *ip;
          continue;
        }
      }

      return oi;
    }
  };

  /*! obtains an `Intersect_2` function object.
   */
  Intersect_2 intersect_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Intersect_2(derived->traits(), derived->intersect_op());
  }
};

//@}

/// `Are_mergeable_2`
//@{

//! Fallback selected `Are_mergeable_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_are_mergeable_2 {};

//! Partial specialization selected if `BaseTraits::Are_mergeable_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_are_mergeable_2<BaseTraits, Derived, std::enable_if_t<has_are_mergeable_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that tests whether two \f$x\f$-monotone curves can be merged.
  class Are_mergeable_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Are_mergeable_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Are_mergeable_2(const Base& base, bool enabled = true) :
      m_object(base.are_mergeable_2_object()), m_enabled(enabled) {}

    /*! determines whether two \f$x\f$-monotone curves can be merged.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the two curve are mergeable and false otherwise.
     * Two curves are mergeable if they have the same underlying theoretical curve.
     */
    bool operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2) const {
      if (! m_enabled) return m_object(xcv1, xcv2);
      std::cout << "are_mergeable" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      bool mergeable = m_object(xcv1, xcv2);
      std::cout << "  result: " << mergeable << std::endl;
      return mergeable;
    }
  };

  /*! obtains an `Are_mergeable_2` function object.
   */
  Are_mergeable_2 are_mergeable_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Are_mergeable_2(derived->traits(), derived->are_mergeable_op());
  }
};

//@}

/// `Merge_2`
//@{

//! Fallback selected `Merge_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_merge_2 {};

//! Partial specialization selected if `BaseTraits::Merge_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_merge_2<BaseTraits, Derived, std::enable_if_t<has_merge_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that merges two \f$x\f$-monotone curves into one.
  class Merge_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Merge_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Merge_2(const Base& base, bool enabled = true) : m_object(base.merge_2_object()), m_enabled(enabled) {}

    /*! merges two \f$x\f$-monotone curves into one.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param xcv the merged curve.
     */
    void operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2, X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv1, xcv2, xcv);
      std::cout << "merge" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      return m_object(xcv1, xcv2, xcv);
      std::cout << "  result: " << xcv << std::endl;
    }
  };

  /*! obtains a `Merge_2` function object.
   */
  Merge_2 merge_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Merge_2(derived->traits(), derived->merge_op());
  }
};

//@}

/// `Construct_opposite_2`
//@{

//! Fallback selected if `Construct_opposite_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_construct_opposite_2 {};

//! Partial specialization selected if `BaseTraits::Construct_opposite_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_construct_opposite_2<BaseTraits, Derived, std::enable_if_t<has_construct_opposite_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A fnuctor that constructs an opposite \f$x\f$-monotone curve.
  class Construct_opposite_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Construct_opposite_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Construct_opposite_2(const Base& base, bool enabled = true) :
      m_object(base.construct_opposite_2_object()), m_enabled(enabled) {}

    /*! constructs an opposite \f$x\f$-monotone curve.
     * \param xcv the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xcv) {
      if (! m_enabled) return m_object(xcv);
      std::cout << "construct_opposite" << std::endl
                << "  xcv: " << xcv << std::endl;
      X_monotone_curve_2 xcv_out = m_object(xcv);
      std::cout << "  result: " << xcv_out << std::endl;
      return xcv;
    }
  };

  /*! obtains a `Construct_opposite_2` function object.
   */
  Construct_opposite_2 construct_opposite_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_opposite_2(derived->traits(), derived->construct_opposite_op());
  }
};

/// `Construct_point_2`
//@{

/* Fallback selected if the functor `BaseTraits::Construct_point_2` does not define an operator that accepts the
 * parameters `const FT&` and `const FT&`.
 */
template <typename, typename, typename = void>
class Tracing_construct_point_2_xy {
protected:
  template <typename T>
  class Construct_point_2 {
  public:
    void operator()() {}; // avoids compilation errors
  };
};

/* Partial specialization selected if the functor `BaseTraits::Construct_point_2` defines an operator that accepts the
 * parameters `const FT&` and `const FT&`.
 */
template <typename BaseTraits, typename Derived>
class Tracing_construct_point_2_xy<BaseTraits, Derived,
                                   std::enable_if_t<has_construct_point_2_xy<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  //! A functor that constructs a point.
  template <typename T>
  class Construct_point_2 {
    using Point_2 = typename Base::Point_2;

  public:
    /*! constructs a point given two coordinates.
     */
    template <typename Coord1, typename Coord2>
    Point_2 operator()(const Coord1& x, const Coord2& y) {
      const T* derived = static_cast<const T*>(this);
      if (! derived->m_enabled) return derived->m_object(x, y);
      std::cout << "construct_point_2" << std::endl
                << "  x: " << x << ", y: " << y << std::endl;
      Point_2 p = derived->m_object(x, y);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };
};

// Fallback selected if `BaseTraits::Construct_point_2` is not defined.
template <typename, typename, typename = void> class Tracing_construct_point_2 {};

// Partial specialization selected if `BaseTraits::Construct_point_2` is defined .
template <typename BaseTraits, typename Derived>
class Tracing_construct_point_2<BaseTraits, Derived, std::enable_if_t<has_construct_point_2<BaseTraits>::value>> :
    public Tracing_construct_point_2_xy<BaseTraits, Derived> {
  using Base = BaseTraits;

public:
  class Construct_point_2; // forward declaration

private:
  using Tracing_construct_point_2_xy =
    typename Tracing_construct_point_2_xy<Base, Derived>::template Construct_point_2<Construct_point_2>;

public:
  //! A functor that constructs a point.
  class Construct_point_2 : Tracing_construct_point_2_xy {
    using Point_2 = typename Base::Point_2;

  private:
    typename Base::Construct_point_2 m_object;
    bool m_enabled;

  public:
    friend Tracing_construct_point_2_xy;

    using Tracing_construct_point_2_xy::operator();

    /*! constructs
     */
    Construct_point_2(const Base& base, bool enabled = true) :
      m_object(base.construct_point_2_object()), m_enabled(enabled) {}

    /*! constructs a point.
     * \return the constructed point.
     */
    template <typename... Args>
    Point_2 operator()(Args... args) const {
      if (! m_enabled) return m_object(std::forward<Args>(args)...);
      std::cout << "construct_point_2" << std::endl
                << "  parameters..." << std::endl;
      Point_2 p = m_object(std::forward<Args>(args)...);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  /*! obtains a `Construct_point_2` function object.
   */
  Construct_point_2 construct_point_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_point_2(derived->traits(), derived->construct_point_op());
  }
};

//@}

/// `Construct_x_monotone_curve_2`
//@{

//! Fallback selected if `Construct_x_monotone_curve_2` is not defined in the base traits
template <typename, typename, typename = void> class Tracing_construct_x_monotone_curve_2 {};

//! Partial specialization selected if `BaseTraits::Construct_x_monotone_curve_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_construct_x_monotone_curve_2<BaseTraits, Derived,
                                           std::enable_if_t<has_construct_x_monotone_curve_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that constructs an \f$x\f$-monotone curve.
  class Construct_x_monotone_curve_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Construct_x_monotone_curve_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Construct_x_monotone_curve_2(const Base& base, bool enabled = true) :
      m_object(base.construct_x_monotone_curve_2_object()), m_enabled(enabled) {}

    /*! constructs an \f$x\f$-monotone curve.
     * \return the constructed \f$x\f$.monotone curve
     */
    template <typename... Args>
    X_monotone_curve_2 operator()(Args... args) const {
      if (! m_enabled) return m_object(std::forward<Args>(args)...);
      std::cout << "construct_x_monotone_curve_2" << std::endl
                << "  parameters..." << std::endl;
      X_monotone_curve_2 xcv = m_object(std::forward<Args>(args)...);
      std::cout << "  result: " << xcv << std::endl;
      return xcv;
    }
  };

  /*! obtains a `Construct_x_monotone_curve_2` function object.
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_x_monotone_curve_2(derived->traits(), derived->construct_point_op());
  }
};

//@}

/// `Construct_curve_2`
//@{

// Fallback selected if `BaseTraits::Construct_curve_2` is not defined.
template <typename, typename, typename = void> class Tracing_construct_curve_2 {};

// Partial specialization selected if `BaseTraits::Construct_curve_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_construct_curve_2<BaseTraits, Derived, std::enable_if_t<has_construct_curve_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that constructs a curve.
  class Construct_curve_2 {
    using Curve_2 = typename Base::Curve_2;

  private:
    typename Base::Construct_curve_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Construct_curve_2(const Base& base, bool enabled = true) :
      m_object(base.construct_curve_2_object()), m_enabled(enabled) {}

    /*! constructs a curve.
     * \return the constructed curve.
     */
    template <typename... Args>
    Curve_2 operator()(Args... args) const {
      if (! m_enabled) return m_object(std::forward<Args>(args)...);
      std::cout << "construct_curve_2" << std::endl
                << "  parameters..." << std::endl;
      Curve_2 cv = m_object(std::forward<Args>(args)...);
      std::cout << "  result: " << cv << std::endl;
      return cv;
    }
  };

  /*! obtains a `Construct_curve_2` function object.
   */
  Construct_curve_2 construct_curve_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_curve_2(derived->traits(), derived->construct_point_op());
  }
};

//@}

/// `Compare_endpoints_xy_2`
//@{

//! Fallback selected if Compare_endpoints_xy_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_compare_endpoints_xy_2 {};

//! Partial specialization selected if `BaseTraits::Compare_endpoints_xy_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_compare_endpoints_xy_2<BaseTraits, Derived,
                                     std::enable_if_t<has_compare_endpoints_xy_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the two endpoints of an \f$x\f$-monotone curve lexigoraphically.
  class Compare_endpoints_xy_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_endpoints_xy_2(const Base& base, bool enabled = true) :
      m_object(base.compare_endpoints_xy_2_object()), m_enabled(enabled) {}

    /*! compares the two endpoints of an \f$x\f$-monotone curve lexigoraphically.
     * \param xcv the curve.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv) {
      if (! m_enabled) return m_object(xcv);
      std::cout << "compare_endpoints_xy" << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains a `Compare_endpoints_xy_2` function object.
   */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_endpoints_xy_2(derived->traits(), derived->compare_endpoints_xy_op());
  }
};

//@}

/// `Approximate_2`
//@{

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts a parameter of type `const Point_2&`.
 */
template <typename, typename, typename = void>
class Tracing_approximate_2_point {
protected:
  template <typename T>
  class Approximate_2 {
  public:
    void operator()() const {}; // avoids compilation errors
  };
};

/* Partial specialization selected if the functor `BaseTraits::Approximate_2` defines an operator that accepts a
 * parameter of type `const Point_2&`.
 */
template <typename BaseTraits, typename Derived>
class Tracing_approximate_2_point<BaseTraits, Derived, std::enable_if_t<has_approximate_2_point<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  using Approximate_point_2 = typename Base::Approximate_point_2;

protected:
  //! A functor that approximates coordinates, points, and \f$x\f$-monotone curves.
  template <typename T>
  class Approximate_2 {
    using Point_2 = typename Base::Point_2;

  public:
    /*! obtains an approximation of a point.
     */
    typename Base::Approximate_point_2 operator()(const Point_2& p) {
      const T* derived = static_cast<const T*>(this);
      if (! derived->m_enabled) return derived->m_object(p);
      std::cout << "approximate" << std::endl
                << "  p: " << p << std::endl;
      auto res = derived->m_object(p);
      std::cout << "  result: " << res << std::endl;
      return res;
    }
  };
};

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts three parameters of the types `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and one optional parameter
 * of type `bool`.
 */
template <typename, typename, typename = void>
class Tracing_approximate_2_xcv {
protected:
  template <typename T>
  class Approximate_2 {
  public:
    void operator()() const {}; // avoids compilation errors
  };
};

/*! Partial specialization selected if the functor `BaseTraits::Approximate_2`
 * defines an operator that accepts three parameters of type `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and one optional parameter
 * of type `bool`.
 */
template <typename BaseTraits, typename Derived>
class Tracing_approximate_2_xcv<BaseTraits, Derived, std::enable_if_t<has_approximate_2_xcv<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  //! A functor that approximates \f$x\f$-monotone curves.
  template <typename T>
  class Approximate_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  public:
    /*! obtains an approximation of an \f$x\f$-monotone curve. */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi, bool l2r = true) {
      const T* derived = static_cast<const T*>(this);
      if (! derived->m_enabled) return derived->m_object(xcv, error, oi, l2r);
      std::cout << "approximate" << std::endl
                << "  xcv: " << xcv << ", error: " << error
                << ", l2r: " << l2r << std::endl;
      std::list<typename Base::Approximate_point_2> container;
      derived->m_object(xcv, error, std::back_inserter(container), l2r);
      if (container.empty()) return oi;

      std::size_t i = 0;
      for (const auto& point : container) {
        std::cout << "  result[" << i++ << "]: " << point << std::endl;
        *oi++ = point;
      }
      return oi;
    }
  };
};

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts four parameters of type `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and `const Bbox_2&`, and
 * one optional parameter of type `bool`.
 */
template <typename, typename, typename = void>
class Tracing_approximate_2_xcv_within_bounds {
protected:
  template <typename T>
  class Approximate_2 {
  public:
    void operator()() const {}; // avoids compilation errors
  };
};

/*! Partial specialization selected if the functor `BaseTraits::Approximate_2`
 * defines an operator that accepts four parameters of type `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and `const Bbox_2&`, and
 * one optional parameter of type `bool`.
 */
template <typename BaseTraits, typename Derived>
class Tracing_approximate_2_xcv_within_bounds<BaseTraits, Derived,
                                              std::enable_if_t<has_approximate_2_xcv_bounds<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  template <typename T>
  class Approximate_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  public:
    /*! obtains an approximation of an \f$x\f$-monotone curve within a given bounding box. */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi, const Bbox_2& bbox,
                              bool l2r = true) const {
      const T* derived = static_cast<const T*>(this);
      if (! derived->m_enabled) return derived->m_object(xcv, error, oi, l2r);
      std::cout << "approximate" << std::endl
                << "  xcv: " << xcv << ", error: " << error
                << ", bbox: " << bbox << ", l2r: " << l2r << std::endl;
      std::list<typename Base::Approximate_point_2> container;
      derived->m_object(xcv, error, std::back_inserter(container), bbox, l2r);
      if (container.empty()) return oi;

      std::size_t i = 0;
      for (const auto& point : container) {
        std::cout << "  result[" << i++ << "]: " << point << std::endl;
        *oi++ = point;
      }
      return oi;
    }
  };
};

// Fallback selected if `Approximate_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_approximate_2 {};

// Partial specialization selected if `BaseTraits::Approximate_2` is defined .
template <typename BaseTraits, typename Derived>
class Tracing_approximate_2<BaseTraits, Derived, std::enable_if_t<has_approximate_2<BaseTraits>::value>> :
    public Tracing_approximate_2_point<BaseTraits, Derived>,
    public Tracing_approximate_2_xcv<BaseTraits, Derived>,
    public Tracing_approximate_2_xcv_within_bounds<BaseTraits, Derived> {
  using Base = BaseTraits;

public:
  class Approximate_2; // forward declaration

private:
  using Tracing_approx_point =
    typename Tracing_approximate_2_point<Base, Derived>::template Approximate_2<Approximate_2>;
  using Tracing_approx_xcv =
    typename Tracing_approximate_2_xcv<Base, Derived>::template Approximate_2<Approximate_2>;
  using Tracing_approx_xcv_within_bounds =
    typename Tracing_approximate_2_xcv_within_bounds<Base, Derived>::template Approximate_2<Approximate_2>;

public:
  using Approximate_number_type = typename Base::Approximate_number_type;

  //! A functor that approximates a coordinates, a point, or an \f$x\f$-monotone curve.
  class Approximate_2 : public Tracing_approx_point,
                        public Tracing_approx_xcv,
                        public Tracing_approx_xcv_within_bounds {
    using Point_2 = typename Base::Point_2;

  public:
    friend Tracing_approx_point;
    friend Tracing_approx_xcv;
    friend Tracing_approx_xcv_within_bounds;

    using Tracing_approx_point::operator();
    using Tracing_approx_xcv::operator();
    using Tracing_approx_xcv_within_bounds::operator();

    /*! constructs
     */
    Approximate_2(const Base& base, bool enabled = true) : m_enabled(enabled), m_object(base.approximate_2_object()) {}

    /*! obtains an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre `i` is either 0 or 1.
     * \return An approximation of `p`'s \f$x\f$-coordinate (if `i` == 0), or an
     *         approximation of `p`'s \f$y\f$-coordinate (if `i` == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) {
      if (! m_enabled) return m_object(p, i);
      std::cout << "approximate" << std::endl
                << "  p: " << p << ", i: " << i << std::endl;
      auto res = m_object(p, i);
      std::cout << "  result: " << res << std::endl;
      return res;
    }

  private:
    bool m_enabled;
    typename Base::Approximate_2 m_object;
  };

  /*! obtains an `Approximate_2` function object.
   */
  Approximate_2 approximate_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Approximate_2(derived->traits(), derived->approximate_op());
  }
};

//@}

/// `Is_on_x_identification_2`
//@{

//! Fallback selected if `Is_on_x_identification_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_is_on_x_identification_2 {};

//! Partial specialization selected if `BaseTraits::Is_on_x_identification_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_is_on_x_identification_2<BaseTraits, Derived,
                                       std::enable_if_t<has_is_on_x_identification_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether a point or curve is on the \f$x\f$-identification curve.
  class Is_on_x_identification_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Is_on_x_identification_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Is_on_x_identification_2(const Base& base, bool enabled = true) :
      m_object(base.is_on_x_identification_2_object()), m_enabled(enabled) {}

    /*! determines whether a point is on the \f$x\f$-identification curve.
     * \param p the point.
     */
    bool operator()(const Point_2& p) const {
      if (! m_enabled) return m_object(p);
      std::cout << "is_on_x_identification" << std::endl
                << "  p: " << p << std::endl;
      bool cr = m_object(p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! determines whether a curve is on the \f$x\f$-identification curve.
     * \param xcv the curve.
     */
    bool operator()(const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv);
      std::cout << "is_on_x_identification" << std::endl
                << "  xcv: " << xcv << std::endl;
      bool cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains an `Is_on_x_identification_2` function object.
   */
  Is_on_x_identification_2 is_on_x_identification_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Is_on_x_identification_2(derived->traits(), derived->is_on_x_identification_op());
  }
};

//@}

/// `is_on_y_identification_2`
//@{

//! Fallback selected if `Is_on_y_identification_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_is_on_y_identification_2 {};

//! Partial specialization selected if `BaseTraits::Is_on_y_identification_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_is_on_y_identification_2<BaseTraits, Derived,
                                       std::enable_if_t<has_is_on_y_identification_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether a point or curve is on the \f$y\f$-identification curve.
  class Is_on_y_identification_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Is_on_y_identification_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Is_on_y_identification_2(const Base& base, bool enabled = true) :
      m_object(base.is_on_y_identification_2_object()), m_enabled(enabled) {}

    /*! determines whether a point s on the \f$y\f$-identification curve.
     * \param p the point.
     */
    bool operator()(const Point_2& p) const {
      if (! m_enabled) return m_object(p);
      std::cout << "is_on_y_identification" << std::endl
                << "  p: " << p << std::endl;
      bool cr = m_object(p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! determines whether a curve is on the \f$y\f$-identification curve.
     * \param xcv the curve.
     */
    bool operator()(const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv);
      std::cout << "is_on_y_identification" << std::endl
                << "  xcv: " << xcv << std::endl;
      bool cr = m_object(xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains an `Is_on_y_identification_2` function object.
   */
  Is_on_y_identification_2 is_on_y_identification_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Is_on_y_identification_2(derived->traits(), derived->is_on_y_identification_op());
  }
};

//@}

/// `Compare_x_on_boundary_2`
//@{

//! Fallback selected if `Compare_x_on_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_compare_x_on_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_x_on_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_compare_x_on_boundary_2<BaseTraits, Derived,
                                      std::enable_if_t<has_compare_x_on_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  /*! A functor that compares the \f$x\f$-coordinate of two given points or
   * curve ends that lie on horizontal boundaries.
   */
  class Compare_x_on_boundary_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_x_on_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_x_on_boundary_2(const Base& base, bool enabled = true) :
      m_object(base.compare_x_on_boundary_2_object()), m_enabled(enabled) {}

    /*! compares the \f$x\f$-coordinate of two given points that lie on horizontal boundaries.
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      if (! m_enabled) return m_object(p1, p2);
      std::cout << "compare_x_on_boundary" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! compares the \f$x\f$-coordinate of a given point and a given curve end that lie on horizontal boundaries.
     * \param pt the point.
     * \param xcv the curve.
     * \param ce the curve-end.
     */
    Comparison_result operator()(const Point_2& pt, const X_monotone_curve_2& xcv, Arr_curve_end ce) {
      if (! m_enabled) return m_object(pt, xcv, ce);
      std::cout << "compare_x_on_boundary" << std::endl
                << "  pt: " << pt << std::endl
                << " xcv: " << xcv << std::endl
                << "  ce: " << ce << std::endl;
      Comparison_result cr = m_object(pt, xcv, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }

    /*! compares the \f$x\f$-coordinate of two given curve ends that lie on horizontal boundaries.
     * \param xcv1 the first curve.
     * \param ce1 the first curve-end.
     * \param xcv2 the second curve.
     * \param ce2 the second curve-end.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2, Arr_curve_end ce2) {
      if (! m_enabled) return m_object(xcv2, ce1, xcv2, ce2);
      std::cout << "compare_x_on_boundary" << std::endl
                << "xcv1: " << xcv1 << std::endl
                << " ce1: " << ce1 << std::endl
                << "xcv2: " << xcv2 << std::endl
                << " ce2: " << ce2 << std::endl;
      Comparison_result cr = m_object(xcv1, ce1, xcv2, ce2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains a `Compare_x_on_boundary_2` function object.
   */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_x_on_boundary_2(derived->traits(), derived->compare_x_on_boundary_op());
  }
};

//@}

/// `Compare_y_on_boundary_2`
//@{

//! Fallback selected if `Compare_y_on_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_compare_y_on_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_y_on_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_compare_y_on_boundary_2<BaseTraits, Derived,
                                      std::enable_if_t<has_compare_y_on_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$y\f$-coordinate of two given points that lie on vertical boundaries.
  class Compare_y_on_boundary_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_y_on_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_y_on_boundary_2(const Base& base, bool enabled = true) :
      m_object(base.compare_y_on_boundary_2_object()),
      m_enabled(enabled)
    {}

    /*! compares the \f$y\f$-coordinate of two given points that lie on vertical boundaries.
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      if (! m_enabled) return m_object(p1, p2);
      std::cout << "compare_y_on_boundary" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains a `Compare_y_on_boundary_2` function object.
   */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_y_on_boundary_2(derived->traits(), derived->compare_y_on_boundary_op());
  }
};

//@}

/// `Compare_x_near_boundary_2`
//@{

//! Fallback selected if `Compare_x_near_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_compare_x_near_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_x_near_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_compare_x_near_boundary_2<BaseTraits, Derived,
                                        std::enable_if_t<has_compare_x_near_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$x\f$-coordinates of curve ends near the boundary of the parameter space.
  class Compare_x_near_boundary_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_x_near_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_x_near_boundary_2(const Base& base, bool enabled = true) :
      m_object(base.compare_x_near_boundary_2_object()), m_enabled(enabled) {}

    /*! compares the \f$x\f$-coordinates of curve ends near the boundary of the parameter space.
     * \param xcv1 the first curve the end of which is to be compared.
     * \param ce1 the identifier of the end of the first curve.
     * \param xcv2 the second curve the end of which is to be compared.
     * \param ce2 the identifier of the end of the second curve.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      if (! m_enabled) return m_object(xcv1, xcv2, ce);
      std::cout << "compare_x_near_boundary" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl
                << "    ce: " << ce << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains a `Compare_x_near_boundary_2` function object.
   */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_x_near_boundary_2(derived->traits(), derived->compare_x_near_boundary_op());
  }
};

//@}

/// `Compare_y_near_boundary_2`
//@{

//! Fallback selected if `Compare_y_near_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Tracing_compare_y_near_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_y_near_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Tracing_compare_y_near_boundary_2<BaseTraits, Derived,
                                        std::enable_if_t<has_compare_y_near_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$y\f$-coordinates of curve ends near the boundary of the parameter space.
  class Compare_y_near_boundary_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_y_near_boundary_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_y_near_boundary_2(const Base& base, bool enabled = true) :
      m_object(base.compare_y_near_boundary_2_object()), m_enabled(enabled) {}

    /*! compares the \f$y\f$-coordinates of curve ends near the boundary of the parameter space.
     * \param xcv1 the first curve the end point of which is tested.
     * \param xcv2 the second curve the end point of which is tested.
     * \param ce the curve-end identifier.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      if (! m_enabled) return m_object(xcv1, xcv2, ce);
      std::cout << "compare_y_near_boundary" << std::endl
                << "  ce: " << ce << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, ce);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /*! obtains a `Compare_y_near_boundary_2` function object.
   */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_y_near_boundary_2(derived->traits(), derived->compare_y_near_boundary_op());
  }
};

} // namespace internal
} // namespace aos2

/*! \class
 * A static metadata traits-class decorator for the arrangement package. It traces the
 * invocations of traits-class functors. It is parameterized with another traits
 * class and inherits from it. For each traits method it prints out its input
 * parameters and its output result.
 *
 * It models all the concepts that the original traits models.
 */
template <typename BaseTraits>
class Arr_tracing_traits_2 :
    public aos2::internal::Tracing_parameter_space_in_x_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_parameter_space_in_y_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_make_x_monotone_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_split_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_do_intersect_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_intersect_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_are_mergeable_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_merge_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_construct_opposite_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_construct_point_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_construct_x_monotone_curve_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_construct_curve_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_compare_endpoints_xy_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_approximate_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_is_on_x_identification_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_is_on_y_identification_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_compare_x_on_boundary_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_compare_y_on_boundary_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_compare_x_near_boundary_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>>,
    public aos2::internal::Tracing_compare_y_near_boundary_2<BaseTraits, Arr_tracing_traits_2<BaseTraits>> {
public:
  enum Operation_id {
    COMPARE_X_2_OP = 0,
    COMPARE_XY_2_OP,
    CONSTRUCT_MIN_VERTEX_2_OP,
    CONSTRUCT_MAX_VERTEX_2_OP,
    IS_VERTICAL_2_OP,
    COMPARE_Y_AT_X_2_OP,
    EQUAL_POINTS_2_OP,
    EQUAL_CURVES_2_OP,
    COMPARE_Y_AT_X_LEFT_2_OP,
    COMPARE_Y_AT_X_RIGHT_2_OP,

    MAKE_X_MONOTONE_2_OP,
    SPLIT_2_OP,
    DO_INTERSECT_2_OP,
    INTERSECT_2_OP,

    ARE_MERGEABLE_2_OP,
    MERGE_2_OP,

    CONSTRUCT_2_OPPOSITE_2_OP,
    CONSTRUCT_POINT_2_OP,
    CONSTRUCT_POINT_2_XY_OP,
    CONSTRUCT_X_MONOTONE_CURVE_2_OP,
    CONSTRUCT_CURVE_2_OP,
    COMPARE_ENDPOINTS_XY_2_OP,
    APPROXIMATE_2_OP,

    PARAMETER_SPACE_IN_X_2_OP,
    IS_ON_X_IDENTIFICATION_2_OP,
    COMPARE_Y_ON_BOUNDARY_2_OP,
    COMPARE_Y_NEAR_BOUNDARY_2_OP,

    PARAMETER_SPACE_IN_Y_2_OP,
    IS_ON_Y_IDENTIFICATION_2_OP,
    COMPARE_X_ON_BOUNDARY_2_OP,
    COMPARE_X_NEAR_BOUNDARY_2_OP,

    NUMBER_OF_OPERATIONS
  };

  using Base = BaseTraits;
  using Shared_base = std::shared_ptr<Base>;

private:
  //! A set of bits that indicate whether operations should be traced.
  unsigned long long m_flags;

  //! The traitse being traced
  Shared_base m_base_traits;

public:
  bool compare_x_op() const { return (0 != (m_flags & (0x1ull << COMPARE_X_2_OP))); }

  bool compare_xy_op() const { return (0 != (m_flags & (0x1ull << COMPARE_XY_2_OP))); }

  bool construct_min_vertex_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_MIN_VERTEX_2_OP))); }

  bool construct_max_vertex_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_MAX_VERTEX_2_OP))); }

  bool is_vertical_op() const { return (0 != (m_flags & (0x1ull << IS_VERTICAL_2_OP))); }

  bool compare_y_at_x_op() const { return (0 != (m_flags & (0x1ull << COMPARE_Y_AT_X_2_OP))); }

  bool equal_points_op() const { return (0 != (m_flags & (0x1ull << EQUAL_POINTS_2_OP))); }

  bool equal_curves_op() const { return (0 != (m_flags & (0x1ull << EQUAL_CURVES_2_OP))); }

  bool compare_y_at_x_left_op() const { return (0 != (m_flags & (0x1ull << COMPARE_Y_AT_X_LEFT_2_OP))); }

  bool compare_y_at_x_right_op() const { return (0 != (m_flags & (0x1ull << COMPARE_Y_AT_X_RIGHT_2_OP))); }

  bool make_x_monotone_op() const { return (0 != (m_flags & (0x1ull << MAKE_X_MONOTONE_2_OP))); }

  bool split_op() const { return (0 != (m_flags & (0x1ull << SPLIT_2_OP))); }

  bool do_intersect_op() const { return (0 != (m_flags & (0x1ull << DO_INTERSECT_2_OP))); }

  bool intersect_op() const { return (0 != (m_flags & (0x1ull << INTERSECT_2_OP))); }

  bool are_mergeable_op() const { return (0 != (m_flags & (0x1ull << ARE_MERGEABLE_2_OP))); }

  bool merge_op() const { return (0 != (m_flags & (0x1ull << MERGE_2_OP))); }

  bool construct_opposite_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_2_OPPOSITE_2_OP))); }

  bool construct_point_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_POINT_2_OP))); }

  bool construct_x_monotone_curve_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_X_MONOTONE_CURVE_2_OP))); }

  bool construct_curve_op() const { return (0 != (m_flags & (0x1ull << CONSTRUCT_CURVE_2_OP))); }

  bool compare_endpoints_xy_op() const { return (0 != (m_flags & (0x1ull << COMPARE_ENDPOINTS_XY_2_OP))); }

  bool approximate_op() const { return (0 != (m_flags & (0x1ull << APPROXIMATE_2_OP))); }

  // left-right

  bool parameter_space_in_x_op() const { return (0 != (m_flags & (0x1ull << PARAMETER_SPACE_IN_X_2_OP))); }

  bool is_on_x_identification_op() const { return (0 != (m_flags & (0x1ull << IS_ON_X_IDENTIFICATION_2_OP))); }

  bool compare_y_on_boundary_op() const { return (0 != (m_flags & (0x1ull << COMPARE_Y_ON_BOUNDARY_2_OP))); }

  bool compare_y_near_boundary_op() const { return (0 != (m_flags & (0x1ull << COMPARE_Y_NEAR_BOUNDARY_2_OP))); }

  // bottom-top

  bool parameter_space_in_y_op() const { return (0 != (m_flags & (0x1ull << PARAMETER_SPACE_IN_Y_2_OP))); }

  bool is_on_y_identification_op() const { return (0 != (m_flags & (0x1ull << IS_ON_Y_IDENTIFICATION_2_OP))); }

  bool compare_x_on_boundary_op() const { return (0 != (m_flags & (0x1ull << COMPARE_X_ON_BOUNDARY_2_OP))); }

  bool compare_x_near_boundary_op() const { return (0 != (m_flags & (0x1ull << COMPARE_X_NEAR_BOUNDARY_2_OP))); }

  /*! constructs default.
   */
  template<typename ... Args>
  Arr_tracing_traits_2(Args ... args) : m_base_traits(std::make_shared<Base>(std::forward<Args>(args)...))
  { enable_all_traces(); }

  /*! constructs from a shared pointer.
   * \param[in] traits the taits being traced.
   * We use std::move to save an atomic operation.  Observe that moving a
   * shared_ptr simply transfers the internal pointer without touching the
   * atomic reference counter. This is faster than a copy.
   */
  Arr_tracing_traits_2(Shared_base traits) : m_base_traits(std::move(traits)) { enable_all_traces(); }

  /*! disables copy constructor.
   */
  Arr_tracing_traits_2(const Arr_tracing_traits_2&) = delete;

  /*! enables the trace of a traits operation.
   * \param id the operation identifier.
   */
  void enable_trace(Operation_id id) { m_flags |= 0x1ull << id; }

  /*! enables the trace of all traits operations.
   */
  void enable_all_traces() { m_flags = 0xffffffff; }

  /*! disables the trace of a traits operation.
   * \param id the operation identifier.
   */
  void disable_trace(Operation_id id) { m_flags &= ~(0x1ull << id); }

  /*! disables the trace of all traits operations.
   */
  void disable_all_traces() { m_flags = 0x0; }

  /*! obtains a const reference to the traits being traced.
   */
  const Base& traits() const { return *m_base_traits; }

  /*! obtains a reference to the traits being traced.
   */
  Base& traits() { return *m_base_traits; }

  /*! obtains the smart pointer to the traits being traced.
   */
  Shared_base shared_traits() const { return m_base_traits; }

  /// \name Types and functors inherited from `BaseTraits`
  //@{

  // Traits types:
  using Has_left_category = typename Base::Has_left_category;
  using Has_merge_category = typename Base::Has_merge_category;

  using Left_side_category = typename internal::Arr_complete_left_side_category<Base>::Category;
  using Bottom_side_category = typename internal::Arr_complete_bottom_side_category<Base>::Category;
  using Top_side_category = typename internal::Arr_complete_top_side_category<Base>::Category;
  using Right_side_category = typename internal::Arr_complete_right_side_category<Base>::Category;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Curve_2 = typename Base::Curve_2;

  //@}

  //! A functor that compares the \f$x\f$-coordinates of two points.
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_x_2(const Base& base, bool enabled = true) :
      m_object(base.compare_x_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param p1 first point.
     * \param p2 second point.
     * \return the comparison result.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      if (! m_enabled) return m_object(p1, p2);
      std::cout << "compare_x" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //! A functor that compares two points lexigoraphically: by \f$x\f$, then by \f$y\f$.
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_xy_2(const Base& base, bool enabled = true) :
      m_object(base.compare_xy_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param p1 the first point.
     * \param p2 the second point.
     * \return the comparison result.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      if (! m_enabled) return m_object(p1, p2);
      std::cout << "compare_xy" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      Comparison_result cr = m_object(p1, p2);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //! A functor that obtains the left endpoint of an \f$x\f$-monotone curve.
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Construct_min_vertex_2(const Base& base, bool enabled = true) :
      m_object(base.construct_min_vertex_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv the curve the left endpoint of which is obtained.
     * \return the left endpoint.
     */
    using Subcurve_ctr_minv = typename Base::Construct_min_vertex_2;
    using Return_type = decltype(std::declval<Subcurve_ctr_minv>().operator()(std::declval<X_monotone_curve_2>()));
    Return_type operator()(const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv);
      std::cout << "construct_min_vertex" << std::endl
                << "  xcv: " << xcv << std::endl;
      Return_type p = m_object(xcv);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  //! A functor that obtains the right endpoint of an \f$x\f$-monotone curve.
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Construct_max_vertex_2(const Base& base, bool enabled = true) :
      m_object(base.construct_max_vertex_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv the curve the right endpoint of which is obtained.
     * \return the right endpoint.
     */
    using Subcurve_ctr_maxv = typename Base::Construct_max_vertex_2;
    using Return_type = decltype(std::declval<Subcurve_ctr_maxv>().operator()(std::declval<X_monotone_curve_2>()));
    Return_type operator()(const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv);
      std::cout << "construct_max_vertex" << std::endl
                << "  xcv: " << xcv << std::endl;
      Return_type p = m_object(xcv);
      std::cout << "  result: " << p << std::endl;
      return p;
    }
  };

  //! A functor that checks whether a given \f$x\f$-monotone curve is vertical.
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Is_vertical_2(const Base& base, bool enabled = true) :
      m_object(base.is_vertical_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv the curve.
     * \return a Boolean that indicates whether the curve is vertical or not.
     */
    bool operator()(const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(xcv);
      std::cout << "is_vertical" << std::endl
                << "  xcv: " << xcv << std::endl;
      bool is_vertical = m_object(xcv);
      std::cout << "  result: " << is_vertical << std::endl;
      return is_vertical;
    }
  };

  /*! A functor that compares the \f$y\f$-coordinates of a point and an
   * \f$x\f$-monotone curve at the point \f$x\f$-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_y_at_x_2(const Base& base, bool enabled = true) :
      m_object(base.compare_y_at_x_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param p the point.
     * \param xcv the curve.
     * \return the comparison result.
     */
    Comparison_result operator()(const Point_2& p,
                                 const X_monotone_curve_2& xcv) const {
      if (! m_enabled) return m_object(p, xcv);
      std::cout << "compare_y_at_x" << std::endl
                << "  p: " << p << std::endl
                << "  xcv: " << xcv << std::endl;
      Comparison_result cr = m_object(p, xcv);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  //! A functor that checks whether two points and two \f$x\f$-monotone curves are identical.
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    bool m_enabled_point;
    bool m_enabled_curve;

  public:
    /*! constructs
     */
    Equal_2(const Base& base, bool enabled_point = true, bool enabled_curve = true) :
      m_object(base.equal_2_object()),
      m_enabled_point(enabled_point),
      m_enabled_curve(enabled_curve)
    {}

    /*! operates
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the \f$x\f$-monotone curves are equal and false
     *         otherwise.
     */
    bool operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2) const {
      if (! m_enabled_curve) return m_object(xcv1, xcv2);
      std::cout << "equal 1" << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv1: " << xcv1 << std::endl;
      bool equal = m_object(xcv1, xcv2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }

    /*! operates
     * \param p1 the first point.
     * \param p2 the second point.
     * \return true if the points are equal and false otherwise.
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      if (! m_enabled_point) return m_object(p1, p2);
      std::cout << "equal 2" << std::endl
                << "  p1: " << p1 << std::endl
                << "  p2: " << p2 << std::endl;
      bool equal = m_object(p1, p2);
      std::cout << "  result: " << equal << std::endl;
      return equal;
    }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the left of their intersection
   * point.
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_y_at_x_left_2(const Base& base, bool enabled = true) :
      m_object(base.compare_y_at_x_left_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param p the reference point.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const {
      if (! m_enabled) return m_object(xcv1, xcv2, p);
      std::cout << "compare_y_at_x_left" << std::endl
                << "  p: " << p << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, p);
      std::cout << "  result:" << cr << std::endl;
      return cr;
    }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the right of their intersection
   * point.
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    bool m_enabled;

  public:
    /*! constructs
     */
    Compare_y_at_x_right_2(const Base& base, bool enabled = true) :
      m_object(base.compare_y_at_x_right_2_object()), m_enabled(enabled) {}

    /*! operates
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param p the reference point.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const {
      if (! m_enabled) return m_object(xcv1, xcv2, p);
      std::cout << "compare_y_at_x_right" << std::endl
                << "  p: " << p << std::endl
                << "  xcv1: " << xcv1 << std::endl
                << "  xcv2: " << xcv2 << std::endl;
      Comparison_result cr = m_object(xcv1, xcv2, p);
      std::cout << "  result: " << cr << std::endl;
      return cr;
    }
  };

  /// \name Obtain the appropriate functor
  //@{

  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(traits(), compare_x_op()); }

  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(traits(), compare_xy_op()); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(traits(), construct_min_vertex_op()); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(traits(), construct_max_vertex_op()); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(traits(), is_vertical_op()); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(traits(), compare_y_at_x_op()); }

  Equal_2 equal_2_object() const
  { return Equal_2(traits(), equal_points_op(), equal_curves_op()); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(traits(), compare_y_at_x_left_op()); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(traits(), compare_y_at_x_right_op()); }

  //@}
};

template <typename OutputStream>
OutputStream& operator<<(OutputStream& os, Comparison_result cr) {
  os << ((cr == SMALLER) ? "SMALLER" : (cr == EQUAL) ? "EQUAL" : "LARGER");
  return os;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
