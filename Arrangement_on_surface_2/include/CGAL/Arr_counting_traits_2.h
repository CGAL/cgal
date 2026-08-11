// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Eric Berberich <ericb@post.tau.ac.il>
//            Shepard Liu    <shepard0liu@gmail.com>

#ifndef CGAL_ARR_COUNTING_TRAITS_H
#define CGAL_ARR_COUNTING_TRAITS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 */

#include <iostream>
#include <atomic>
#include <array>
#include <memory>

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_has.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {
namespace aos2 {
namespace internal {

/// `Parameter_space_in_x_2`
//@{

//! Fallback selected `Parameter_space_in_x_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_parameter_space_in_x_2 {};

//! Partial specialization selected if `BaseTraits::`Parameter_space_in_x_2 is defined.
template <typename BaseTraits, typename Derived>
class Counting_parameter_space_in_x_2<BaseTraits, Derived,
                                      std::enable_if_t<has_parameter_space_in_x_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone curve
  * lies on a boundary of the parameter space along the \f$x\f$-axis.
  */
  class Parameter_space_in_x_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Parameter_space_in_x_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs
     */
    Parameter_space_in_x_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.parameter_space_in_x_2_object()), m_counter1(counter1), m_counter2(counter2) {}

    /*! operates
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc, Arr_curve_end ce) const {
      ++m_counter1;
      return m_object(xc, ce);
    }

    /*! operates
     */
    Arr_parameter_space operator()(const Point_2& p) const {
      ++m_counter2;
      return m_object(p);
    }
  };

  /*! obtains a `Parameter_space_in_x_2` function object.
   */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Parameter_space_in_x_2(derived->traits(),
                                  derived->m_counters[Derived::PARAMETER_SPACE_IN_X_2_CURVE_END_OP],
                                  derived->m_counters[Derived::PARAMETER_SPACE_IN_X_2_POINT_OP]);
  }
};

//@}

/// `Parameter_space_in_y_2`
//@{

//! Fallback selected `Parameter_space_in_y_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_parameter_space_in_y_2 {};

//! Partial specialization selected if `BaseTraits::Parameter_space_in_y_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_parameter_space_in_y_2<BaseTraits, Derived,
                                      std::enable_if_t<has_parameter_space_in_y_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  /*! A functor that determines whether an endpoint of an \f$x\f$-monotone arc
   * lies on a boundary of the parameter space along the \f$y\f$-axis.
   */
  class Parameter_space_in_y_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Parameter_space_in_y_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs
     */
    Parameter_space_in_y_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.parameter_space_in_y_2_object()), m_counter1(counter1), m_counter2(counter2) {}

    /*! operates
     */
    Arr_parameter_space operator()(const X_monotone_curve_2& xc, Arr_curve_end ce) const {
      ++m_counter1;
      return m_object(xc, ce);
    }

    /*! operates
     */
    Arr_parameter_space operator()(const Point_2& p) const {
      ++m_counter2;
      return m_object(p);
    }
  };

  /*! obtains a `Parameter_space_in_y_2` function object.
   */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Parameter_space_in_y_2(derived->traits(),
                                  derived->m_counters[Derived::PARAMETER_SPACE_IN_Y_2_CURVE_END_OP],
                                  derived->m_counters[Derived::PARAMETER_SPACE_IN_Y_2_POINT_OP]);
  };
};

//@}

/// `Make_x_monotone_2`
//@{

//! Fallback selected `Make_x_monotone_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_make_x_monotone_2 {};

//! Partial specialization selected if `BaseTraits::Make_x_monotone_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_make_x_monotone_2<BaseTraits, Derived, std::enable_if_t<has_make_x_monotone_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  using Curve_2 = typename Base::Curve_2;

  //! A functor that subdivides a curve into \f$x\f$-monotone curves.
  class Make_x_monotone_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Make_x_monotone_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Make_x_monotone_2(const Base& base, std::size_t& counter) :
      m_object(base.make_x_monotone_2_object()), m_counter(counter) {}

    /*! subdivides a given curve into \f$x\f$-monotone subcurves and insert them into a given output iterator.
     * \param cv the curve.
     * \param oi the output iterator for the result. Its value type is a variant
     *           that wraps `Point_2` or an `X_monotone_curve_2` objects.
     * \return The past-the-end iterator.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
      ++m_counter;
      return m_object(cv, oi);
    }
  };

  /*! obtains a `Make_x_monotone_2` function object.
   */
  Make_x_monotone_2 make_x_monotone_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Make_x_monotone_2(derived->traits(), derived->m_counters[Derived::MAKE_X_MONOTONE_2_OP]);
  }
};

//@}

/// `Split_2`
//@{

//! Fallback selected `Split_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_split_2 {};

//! Partial specialization selected if `BaseTraits::Split_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_split_2<BaseTraits, Derived, std::enable_if_t<has_split_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that splits an \f$x\f$-monotone curve at a point.
  class Split_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Split_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Split_2(const Base& base, std::size_t& counter) : m_object(base.split_2_object()), m_counter(counter) {}

    /*! operates
     */
    void operator()(const X_monotone_curve_2& xc, const Point_2& p,
                    X_monotone_curve_2& xcv1, X_monotone_curve_2& xcv2) const {
      ++m_counter;
      m_object(xc, p, xcv1, xcv2);
    }
  };

  /*! obtains a `Split_2` function object.
   */
  Split_2 split_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Split_2(derived->traits(), derived->m_counters[Derived::SPLIT_2_OP]);
  }
};

//@}

/// `Do_intersect_2`
//@{

//! Fallback selected `Do_intersect_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_do_intersect_2 {};

//! Partial specialization selected if `BaseTraits::Do_intersect_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_do_intersect_2<BaseTraits, Derived, std::enable_if_t<has_do_intersect_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether two \f$x\f$-monotone curves intersect.
  class Do_intersect_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
    using Point_2 = typename Base::Point_2;

  private:
    typename Base::Do_intersect_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Do_intersect_2(const Base& base, std::size_t& counter) :
      m_object(base.do_intersect_2_object()), m_counter(counter) {}

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
      ++m_counter;
      return m_object(xcv1, xcv2, consider_common_endpoints);
    }
  };

  /*! obtains a `Do_intersect_2` function object.
   */
  Do_intersect_2 do_intersect_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Do_intersect_2(derived->traits(), derived->m_counters[Derived::DO_INTERSECT_2_OP]);
  }
};

//@}

/// `Intersect_2`
//@{

//! Fallback selected `Intersect_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_intersect_2 {};

//! Partial specialization selected if `BaseTraits::Intersect_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_intersect_2<BaseTraits, Derived, std::enable_if_t<has_intersect_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  using Multiplicity = typename Base::Multiplicity;

  //! A functor that computes intersections between \f$x\f$-monotone curves.
  class Intersect_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
    using Point_2 = typename Base::Point_2;

  private:
    typename Base::Intersect_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Intersect_2(const Base& base, std::size_t& counter) : m_object(base.intersect_2_object()), m_counter(counter) {}

    /*! operates
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2, OutputIterator oi) const {
      ++m_counter;
      return m_object(xcv1, xcv2, oi);
    }
  };

  /*! obtains an `Intersect_2` function object.
   */
  Intersect_2 intersect_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Intersect_2(derived->traits(), derived->m_counters[Derived::INTERSECT_2_OP]);
  }
};

//@}

/// `Are_mergeable_2`
//@{

//! Fallback selected `Are_mergeable_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_are_mergeable_2 {};

//! Partial specialization selected if `BaseTraits::Are_mergeable_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_are_mergeable_2<BaseTraits, Derived, std::enable_if_t<has_are_mergeable_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that tests whether two \f$x\f$-monotone curves can be merged.
  class Are_mergeable_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Are_mergeable_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Are_mergeable_2(const Base& base, std::size_t& counter) :
      m_object(base.are_mergeable_2_object()), m_counter(counter) {}

    /*! determines whether two \f$x\f$-monotone curves can be merged.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \return true if the two curve are mergeable and false otherwise.
     * Two curves are mergeable if they have the same underlying theoretical curve.
     */
    bool operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2) const {
      ++m_counter;
      return m_object(xcv1, xcv2);
    }
  };

  /*! obtains an `Are_mergeable_2` function object.
   */
  Are_mergeable_2 are_mergeable_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Are_mergeable_2(derived->traits(), derived->m_counters[Derived::ARE_MERGEABLE_2_OP]);
  }
};

//@}

/// `Merge_2`
//@{

//! Fallback selected `Merge_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_merge_2 {};

//! Partial specialization selected if `BaseTraits::Merge_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_merge_2<BaseTraits, Derived, std::enable_if_t<has_merge_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that merges two \f$x\f$-monotone curves into one.
  class Merge_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Merge_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Merge_2(const Base& base, std::size_t& counter) : m_object(base.merge_2_object()), m_counter(counter) {}

    /*! merges two \f$x\f$-monotone curves into one.
     * \param xcv1 the first curve.
     * \param xcv2 the second curve.
     * \param xcv the merged curve.
     */
    void operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2, X_monotone_curve_2& xcv) const {
      ++m_counter;
      m_object(xcv1, xcv2, xcv);
    }
  };

  /*! obtains a `Merge_2` function object.
   */
  Merge_2 merge_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Merge_2(derived->traits(), derived->m_counters[Derived::MERGE_2_OP]);
  }
};

//@}

/// `Construct_opposite_2`
//@{

//! Fallback selected if `Construct_opposite_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_construct_opposite_2 {};

//! Partial specialization selected if `BaseTraits::Construct_opposite_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_construct_opposite_2<BaseTraits, Derived,
                                    std::enable_if_t<has_construct_opposite_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that constructs an opposite \f$x\f$-monotone curve.
  class Construct_opposite_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Construct_opposite_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Construct_opposite_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_opposite_2_object()), m_counter(counter) {}

    /*! constructs an opposite \f$x\f$-monotone curve.
     * \param xcv the curve.
     * \return the opposite curve.
     */
    X_monotone_curve_2 operator()(const X_monotone_curve_2& xc) const {
      ++m_counter;
      return m_object(xc);
    }
  };

  /*! obtains a `Construct_opposite_2` function object.
   */
  Construct_opposite_2 construct_opposite_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_opposite_2(derived->traits(), derived->m_counters[Derived::CONSTRUCT_2_OPPOSITE_2_OP]);
  }
};

//@}

/// `Construct_point_2`
//@{

/*! Fallback selected if the functor `BaseTraits::Construct_point_2` does not define an operator that accepts two
 * parameters of type `const FT&`.
 */
template <typename, typename, typename = void>
class Counting_construct_point_2_xy {
protected:
  template <typename T>
  class Construct_point_2 {
  public:
    void operator()() const {}; // avoids compilation errors
  };
};

/*! Partial specialization selected if the functor `BaseTraits::Construct_point_2` defines an operator that accepts two
 * parameters of type `const FT&`.
 */
template <typename BaseTraits, typename Derived>
class Counting_construct_point_2_xy<BaseTraits, Derived,
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
    Point_2 operator()(const Coord1& x, const Coord2& y) const {
      const T* derived = static_cast<const T*>(this);
      ++derived->m_counter2;
      return derived->m_object(x, y);
    }
  };
};

//! Fallback selected `Construct_point_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_construct_point_2 {};

//! Partial specialization selected if `BaseTraits::Construct_point_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_construct_point_2<BaseTraits, Derived, std::enable_if_t<has_construct_point_2<BaseTraits>::value>> :
    public Counting_construct_point_2_xy<BaseTraits, Derived> {
  using Base = BaseTraits;

public:
  class Construct_point_2; // forward declaration

private:
  using Counting_construct_point_2_xy =
    typename Counting_construct_point_2_xy<Base, Derived>::template Construct_point_2<Construct_point_2>;

public:
  //! A functor that constructs a point.
  class Construct_point_2 : public Counting_construct_point_2_xy {
    using Point_2 = typename Base::Point_2;

  private:
    typename Base::Construct_point_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    friend Counting_construct_point_2_xy;

    using Counting_construct_point_2_xy::operator();

    /*! constructs
     */
    Construct_point_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.construct_point_2_object()), m_counter1(counter1), m_counter2(counter2) {}

    /*! constructs a point.
     * \return the constructed point.
     */
    template <typename... Args>
    Point_2 operator()(Args... args) const {
      ++m_counter1;
      return m_object(std::forward<Args>(args)...);
    }
  };

  /*! obtains a `Construct_point_2` function object.
   */
  Construct_point_2 construct_point_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_point_2(derived->traits(),
                             derived->m_counters[Derived::CONSTRUCT_POINT_2_OP],
                             derived->m_counters[Derived::CONSTRUCT_POINT_2_XY_OP]);
  }
};

//@}

/// `Construct_x_monotone_curve_2`
//@{

//! Fallback selected if `Construct_x_monotone_curve_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_construct_x_monotone_curve_2 {};

//! Partial specialization selected if `BaseTraits::Construct_x_monotone_curve_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_construct_x_monotone_curve_2<BaseTraits, Derived,
                                            std::enable_if_t<has_construct_x_monotone_curve_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that constructs an \f$x\f$-monotone curve.
  class Construct_x_monotone_curve_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Construct_x_monotone_curve_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Construct_x_monotone_curve_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_x_monotone_curve_2_object()), m_counter(counter) {}

    /*! constructs an \f$x\f$-monotone curve.
     * \return the constructed \f$x\f$.monotone curve
     */
    template <typename... Args>
    X_monotone_curve_2 operator()(Args... args) const {
      ++m_counter;
      return m_object(std::forward<Args>(args)...);
    }
  };

  /*! obtains a `Construct_x_monotone_curve_2` function object.
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_x_monotone_curve_2(derived->traits(),
                                        derived->m_counters[Derived::CONSTRUCT_X_MONOTONE_CURVE_2_OP]);
  }
};

//@}

/// `Construct_curve_2`
//@{

//! Fallback selected if `Construct_curve_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_construct_curve_2 {};

//! Partial specialization selected if `BaseTraits::Construct_curve_2` is defined
template <typename BaseTraits, typename Derived>
class Counting_construct_curve_2<BaseTraits, Derived, std::enable_if_t<has_construct_curve_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that constructs a curve.
  class Construct_curve_2 {
    /* Defining `Curve_2` in the outer class `Counting_construct_curve_2` would
     * cause an ambiguous definition because it is already defined in the
     * (outer) class `Counting_make_x_monotone_2`. Thus, we only define it in the
     * inner class `Construct_curve_2`. If the base traits class does not
     * support `Make_x_monotone_2`, `Curve_2` will end up being
     * undefined. However, the scenario where the base traits class supports
     * `Construct_curve_2` but does not support `Make_x_monotone_2` is
     * unlikely. If this scenario is encountered after all, a different solution
     * must be devised.
     */
    using Curve_2 = typename Base::Curve_2;

  private:
    typename Base::Construct_curve_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Construct_curve_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_curve_2_object()), m_counter(counter) {}

    /*! constructs a curve.
     * \return the constructed curve.
     */
    template <typename... Args>
    Curve_2 operator()(Args... args) const {
      ++m_counter;
      return m_object(std::forward<Args>(args)...);
    }
  };

  /*! obtains a `Construct_curve_2` function object.
   */
  Construct_curve_2 construct_curve_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Construct_curve_2(derived->traits(), derived->m_counters[Derived::CONSTRUCT_CURVE_2_OP]);
  }
};

//@}

/// `Compare_endpoints_xy_2`
//@{

//! Fallback selected if Compare_endpoints_xy_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_compare_endpoints_xy_2 {};

//! Partial specialization selected if `BaseTraits::Compare_endpoints_xy_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_compare_endpoints_xy_2<BaseTraits, Derived,
                                      std::enable_if_t<has_compare_endpoints_xy_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the two endpoints of an \f$x\f$-monotone curve lexigoraphically.
  class Compare_endpoints_xy_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_endpoints_xy_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_endpoints_xy_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_endpoints_xy_2_object()), m_counter(counter) {}

    /*! compares the two endpoints of an \f$x\f$-monotone curve lexigoraphically.
     * \param xcv the curve.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xc) const {
      ++m_counter;
      return m_object(xc);
    }
  };

  /*! obtains a `Compare_endpoints_xy_2` function object.
   */
  Compare_endpoints_xy_2 compare_endpoints_xy_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_endpoints_xy_2(derived->traits(), derived->m_counters[Derived::COMPARE_ENDPOINTS_XY_2_OP]);
  }
};

//@}

/// \name `Approximate_2`
//@{

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts a parameter of type `const Point_2&`.
 */
template <typename, typename, typename = void>
class Counting_approximate_2_point {
protected:
  template <typename T>
  class Approximate_2 {
  public:
    void operator()() const {}; // avoids compilation errors
  };
};

/*! Partial specialization selected if the functor `BaseTraits::Approximate_2` defines an operator that accepts a
 * parameter of type `const Point_2&`.
 */
template <typename BaseTraits, typename Derived>
class Counting_approximate_2_point<BaseTraits, Derived, std::enable_if_t<has_approximate_2_point<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  //! A functor that approximates coordinates, points, and \f$x\f$-monotone curves.
  template <typename T>
  class Approximate_2 {
    using Point_2 = typename Base::Point_2;
    using Base_approximate_2 = typename Base::Approximate_2;

  public:
    using Approximate_kernel = typename Base_approximate_2::Approximate_kernel;
    using Approximate_point_2 = typename Base_approximate_2::Approximate_point_2;

    /*! obtains an approximation of a point.
     */
    Approximate_point_2 operator()(const Point_2& p) const {
      const T* derived = static_cast<const T*>(this);
      ++derived->m_counter2;
      return derived->m_object(p);
    }
  };
};

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts three parameters of the types `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and one optional parameter
 * of type `bool`.
 */
template <typename, typename, typename = void>
class Counting_approximate_2_xcv {
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
class Counting_approximate_2_xcv<BaseTraits, Derived, std::enable_if_t<has_approximate_2_xcv<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  //! A functor that approximates \f$x\f$-monotone curves.
  template <typename T>
  class Approximate_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  public:
    /*! obtains an approximation of an \f$x\f$-monotone curve. */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi, bool l2r = true) const {
      const T* derived = static_cast<const T*>(this);
      ++derived->m_counter3;
      return derived->m_object(xcv, error, oi, l2r);
    }
  };
};

/*! Fallback selected if the functor `Approximate_2`, nested in the base traits,
 * does not define an operator that accepts four parameters of type `const
 * X_monotone_curve_2&`, `double`, `OutputIterator`, and `const Bbox_2&`, and
 * one optional parameter of type `bool`.
 */
template <typename, typename, typename = void>
class Counting_approximate_2_xcv_within_bounds {
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
class Counting_approximate_2_xcv_within_bounds<BaseTraits, Derived,
                                               std::enable_if_t<has_approximate_2_xcv_bounds<BaseTraits>::value>> {
  using Base = BaseTraits;

protected:
  //! A functor that approximates \f$x\f$-monotone curves within bounds.
  template <typename T>
  class Approximate_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  public:
    /*! obtains an approximation of an \f$x\f$-monotone curve within bounds.
     */
    template <typename OutputIterator>
    OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi, const Bbox_2& bbox,
                              bool l2r = true) const {
      const T* derived = static_cast<const T*>(this);
      ++derived->m_counter4;
      return derived->m_object(xcv, error, oi, bbox, l2r);
    }
  };
};

//! Fallback selected if `Approximate_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_approximate_2 {};

//! Partial specialization selected if `BaseTraits::Approximate_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_approximate_2<BaseTraits, Derived, std::enable_if_t<has_approximate_2<BaseTraits>::value>> :
    public Counting_approximate_2_point<BaseTraits, Derived>,
    public Counting_approximate_2_xcv<BaseTraits, Derived>,
    public Counting_approximate_2_xcv_within_bounds<BaseTraits, Derived> {
  using Base = BaseTraits;

public:
  class Approximate_2; // forward declaration

private:
  using Counting_approx_point =
    typename Counting_approximate_2_point<Base, Derived>::template Approximate_2<Approximate_2>;
  using Counting_approx_xcv =
    typename Counting_approximate_2_xcv<Base, Derived>::template Approximate_2<Approximate_2>;
  using Counting_approx_xcv_within_bounds =
    typename Counting_approximate_2_xcv_within_bounds<Base, Derived>::template Approximate_2<Approximate_2>;

public:
  //! A functor that approximates a coordinates, a point, or an \f$x\f$-monotone curve.
  class Approximate_2 : public Counting_approx_point,
                        public Counting_approx_xcv,
                        public Counting_approx_xcv_within_bounds {
    using Point_2 = typename Base::Point_2;
    using Base_approximate_2 = typename Base::Approximate_2;

  public:
    friend Counting_approx_point;
    friend Counting_approx_xcv;
    friend Counting_approx_xcv_within_bounds;

    using Counting_approx_point::operator();
    using Counting_approx_xcv::operator();
    using Counting_approx_xcv_within_bounds::operator();

    using Approximate_number_type = typename Base_approximate_2::Approximate_number_type;

    /*! constructs
     */
    Approximate_2(const Base& base, std::size_t& counter1, std::size_t& counter2,
                  std::size_t& counter3, std::size_t& counter4) :
      m_object(base.approximate_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3),
      m_counter4(counter4)
    {}

    /*! obtains an approximation of a point coordinate.
     * \param p the exact point.
     * \param i the coordinate index (either 0 or 1).
     * \pre `i` is either 0 or 1.
     * \return An approximation of `p`'s \f$x\f$-coordinate (if `i` == 0), or an
     *         approximation of `p`'s \f$y\f$-coordinate (if `i` == 1).
     */
    Approximate_number_type operator()(const Point_2& p, int i) const {
      ++m_counter1;
      return m_object(p, i);
    }

  private:
    Base_approximate_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;
    std::size_t& m_counter3;
    std::size_t& m_counter4;
  };

  /*! obtains an `Approximate_2` function object.
   */
  Approximate_2 approximate_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Approximate_2(derived->traits(),
                         derived->m_counters[Derived::APPROXIMATE_2_COORD_OP],
                         derived->m_counters[Derived::APPROXIMATE_2_POINT_OP],
                         derived->m_counters[Derived::APPROXIMATE_2_CURVE_OP],
                         derived->m_counters[Derived::APPROXIMATE_2_BOUNDED_CURVE_OP]);
  }
};

//@}

/// `Is_on_x_identification_2`
//@{

//! Fallback selected if `Is_on_x_identification_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_is_on_x_identification_2 {};

//! Partial specialization selected if `BaseTraits::Is_on_x_identification_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_is_on_x_identification_2<BaseTraits, Derived,
                                        std::enable_if_t<has_is_on_x_identification_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether a point or a curve lies on \f$x\f$-identification curve.
  class Is_on_x_identification_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Is_on_x_identification_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs
     */
    Is_on_x_identification_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.is_on_x_identification_2_object()), m_counter1(counter1), m_counter2(counter2) {}

    /*! determines whether a point is on the \f$x\f$-identification curve.
     * \param p the point.
     */
    bool operator()(const Point_2& p) const {
      ++m_counter1;
      return m_object(p);
    }

    /*! determines whether a curve is on the \f$x\f$-identification curve.
     * \param xcv the curve.
     */
    bool operator()(const X_monotone_curve_2& xc) const {
      ++m_counter2;
      return m_object(xc);
    }
  };

  /*! obtains an `Is_on_x_identification_2` function object.
   */
  Is_on_x_identification_2 is_on_x_identification_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Is_on_x_identification_2(derived->traits(),
                                    derived->m_counters[Derived::IS_ON_X_IDENTIFICATION_POINT_2_OP],
                                    derived->m_counters[Derived::IS_ON_X_IDENTIFICATION_CURVE_2_OP]);
  }
};

//@}

/// `is_on_y_identification_2`
//@{

//! Fallback selected if `Is_on_y_identification_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_is_on_y_identification_2 {};

//! Partial specialization selected if `BaseTraits::Is_on_y_identification_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_is_on_y_identification_2<BaseTraits, Derived,
                                        std::enable_if_t<has_is_on_y_identification_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that determines whether a point or a curve lies on the \f$y\f$-identification curve.
  class Is_on_y_identification_2 {
    using Point_2 = typename Base::Point_2;
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Is_on_y_identification_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs
     */
    Is_on_y_identification_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.is_on_y_identification_2_object()), m_counter1(counter1),  m_counter2(counter2) {}

    /*! determines whether a point is on the \f$y\f$-identification curve.
     * \param p the point.
     */
    bool operator()(const Point_2& p) const {
      ++m_counter1;
      return m_object(p);
    }

    /*! determines whether a curve is on the \f$y\f$-identification curve.
     * \param p the point.
     */
    bool operator()(const X_monotone_curve_2& xc) const {
      ++m_counter2;
      return m_object(xc);
    }
  };

  /*! obtains an `Is_on_y_identification_2` function object.
   */
  Is_on_y_identification_2 is_on_y_identification_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Is_on_y_identification_2(derived->traits(),
                                    derived->m_counters[Derived::IS_ON_Y_IDENTIFICATION_2_POINT_OP],
                                    derived->m_counters[Derived::IS_ON_Y_IDENTIFICATION_2_CURVE_OP]);
  }
};

//@}

/// `Compare_x_on_boundary_2`
//@{

//! Fallback selected if `Compare_x_on_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_x_on_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_x_on_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_x_on_boundary_2<BaseTraits, Derived, std::enable_if_t<has_compare_x_on_boundary_2<BaseTraits>::value>> {
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
    std::size_t& m_counter1;
    std::size_t& m_counter2;
    std::size_t& m_counter3;

  public:
    /*! constructs
     */
    Compare_x_on_boundary_2(const Base& base,  std::size_t& counter1, std::size_t& counter2, std::size_t& counter3) :
      m_object(base.compare_x_on_boundary_2_object()),
      m_counter1(counter1),
      m_counter2(counter2),
      m_counter3(counter3)
    {}

    /*! compares the \f$x\f$-coordinate of two given points that lie on horizontal boundaries.
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      ++m_counter1;
      return m_object(p1, p2);
    }

    /*! compares the \f$x\f$-coordinate of a given point and a given curve end that lie on horizontal boundaries.
     * \param pt the point.
     * \param xcv the curve.
     * \param ce the curve-end.
     */
    Comparison_result operator()(const Point_2& pt, const X_monotone_curve_2& xcv, Arr_curve_end ce) const {
      ++m_counter2;
      return m_object(pt, xcv, ce);
    }

    /*! compares the \f$x\f$-coordinate of two given curve ends that lie on horizontal boundaries.
     * \param xcv1 the first curve.
     * \param ce1 the first curve-end.
     * \param xcv2 the second curve.
     * \param ce2 the second curve-end.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, Arr_curve_end ce1,
                                 const X_monotone_curve_2& xcv2, Arr_curve_end ce2) const {
      ++m_counter3;
      return m_object(xcv1, ce1, xcv2, ce2);
    }
  };

  /*! obtains a `Compare_x_on_boundary_2` function object.
   */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_x_on_boundary_2(derived->traits(),
                                   derived->m_counters[Derived::COMPARE_X_ON_BOUNDARY_2_POINTS_OP],
                                   derived->m_counters[Derived::COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP],
                                   derived->m_counters[Derived::COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP]);
  }
};

//@}

/// `Compare_y_on_boundary_2`
//@{

//! Fallback selected if `Compare_y_on_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_compare_y_on_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_y_on_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_compare_y_on_boundary_2<BaseTraits, Derived,
                                       std::enable_if_t<has_compare_y_on_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$y\f$-coordinate of two given points that lie on vertical boundaries.
  class Compare_y_on_boundary_2 {
    using Point_2 = typename Base::Point_2;

  private:
    typename Base::Compare_y_on_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_y_on_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_on_boundary_2_object()), m_counter(counter) {}

    /*! compares the \f$y\f$-coordinate of two given points that lie on vertical boundaries.
     * \param p1 the first point.
     * \param p2 the second point.
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      ++m_counter;
      return m_object(p1, p2);
    }
  };

  /*! obtains a `Compare_y_on_boundary_2` function object.
   */
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_y_on_boundary_2(derived->traits(), derived->m_counters[Derived::COMPARE_Y_ON_BOUNDARY_2_OP]);
  }
};

//@}

/// `Compare_x_near_boundary_2`
//@{

//! Fallback selected if `Compare_x_near_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_compare_x_near_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_x_near_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_compare_x_near_boundary_2<BaseTraits, Derived,
                                         std::enable_if_t<has_compare_x_near_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$x\f$-coordinates of curve ends near the boundary of the parameter space.
  class Compare_x_near_boundary_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_x_near_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_x_near_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_x_near_boundary_2_object()), m_counter(counter) {}

    /*! compares the \f$x\f$-coordinates of curve ends near the boundary of the parameter space.
     * \param xcv1 the first curve the end of which is to be compared.
     * \param ce1 the identifier of the end of the first curve.
     * \param xcv2 the second curve the end of which is to be compared.
     * \param ce2 the identifier of the end of the second curve.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      ++m_counter;
      return m_object(xcv1, xcv2, ce);
    }
  };

  /*! obtains a `Compare_x_near_boundary_2` function object.
   */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_x_near_boundary_2(derived->traits(), derived->m_counters[Derived::COMPARE_X_NEAR_BOUNDARY_2_OP]);
  }
};

//@}

/// `Compare_y_near_boundary_2`
//@{

//! Fallback selected if `Compare_y_near_boundary_2` is not defined in the base traits.
template <typename, typename, typename = void> class Counting_compare_y_near_boundary_2 {};

//! Partial specialization selected if `BaseTraits::Compare_y_near_boundary_2` is defined.
template <typename BaseTraits, typename Derived>
class Counting_compare_y_near_boundary_2<BaseTraits, Derived,
                                         std::enable_if_t<has_compare_y_near_boundary_2<BaseTraits>::value>> {
  using Base = BaseTraits;

public:
  //! A functor that compares the \f$y\f$-coordinates of curve ends near the boundary of the parameter space.
  class Compare_y_near_boundary_2 {
    using X_monotone_curve_2 = typename Base::X_monotone_curve_2;

  private:
    typename Base::Compare_y_near_boundary_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_y_near_boundary_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_near_boundary_2_object()), m_counter(counter) {}

    /*! compares the \f$y\f$-coordinates of curve ends near the boundary of the parameter space.
     * \param xcv1 the first curve the end point of which is tested.
     * \param xcv2 the second curve the end point of which is tested.
     * \param ce the curve-end identifier.
     * \return the comparison result.
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 Arr_curve_end ce) const {
      ++m_counter;
      return m_object(xcv1, xcv2, ce);
    }
  };

  /*! obtains a `Compare_y_near_boundary_2` function object.
   */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const {
    const Derived* derived = static_cast<const Derived*>(this);
    return Compare_y_near_boundary_2(derived->traits(), derived->m_counters[Derived::COMPARE_Y_NEAR_BOUNDARY_2_OP]);
  }
};

//@}

}
}

/*! \class
 * A metadata traits-class decorator for the arrangement package. It counts the
 * number of invocations of traits-class functors. It is parameterized with
 * another traits class and inherits from it. For each traits method it
 * maintains a counter that counts the number of invocations into the method.
 *
 * It models all the concepts that the original traits models.
 */
template <typename BaseTraits>
class Arr_counting_traits_2 :
    public aos2::internal::Counting_parameter_space_in_x_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_parameter_space_in_y_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_make_x_monotone_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_split_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_do_intersect_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_intersect_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_are_mergeable_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_merge_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_construct_opposite_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_construct_point_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_construct_x_monotone_curve_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_construct_curve_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_compare_endpoints_xy_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_approximate_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_is_on_x_identification_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_is_on_y_identification_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_compare_y_on_boundary_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_compare_y_near_boundary_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_x_on_boundary_2<BaseTraits, Arr_counting_traits_2<BaseTraits>>,
    public aos2::internal::Counting_compare_x_near_boundary_2<BaseTraits, Arr_counting_traits_2<BaseTraits>> {

  using Base = BaseTraits;

  using Counting_parameter_space_in_x_2 =
    aos2::internal::Counting_parameter_space_in_x_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_parameter_space_in_y_2 =
    aos2::internal::Counting_parameter_space_in_y_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_make_x_monotone_2 =
    aos2::internal::Counting_make_x_monotone_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_split_2 =
    aos2::internal::Counting_split_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_do_intersect_2 =
    aos2::internal::Counting_do_intersect_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_intersect_2 =
    aos2::internal::Counting_intersect_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_are_mergeable_2 =
    aos2::internal::Counting_are_mergeable_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_merge_2 =
    aos2::internal::Counting_merge_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_construct_opposite_2 =
    aos2::internal::Counting_construct_opposite_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_construct_point_2 =
    aos2::internal::Counting_construct_point_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_construct_x_monotone_curve_2 =
    aos2::internal::Counting_construct_x_monotone_curve_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_construct_curve_2 =
    aos2::internal::Counting_construct_curve_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_compare_endpoints_xy_2 =
    aos2::internal::Counting_compare_endpoints_xy_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_approximate_2 =
    aos2::internal::Counting_approximate_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_is_on_x_identification_2 =
    aos2::internal::Counting_is_on_x_identification_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_is_on_y_identification_2 =
    aos2::internal::Counting_is_on_y_identification_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_compare_y_on_boundary_2 =
    aos2::internal::Counting_compare_y_on_boundary_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_compare_y_near_boundary_2 =
    aos2::internal::Counting_compare_y_near_boundary_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_x_on_boundary_2 =
    aos2::internal::Counting_x_on_boundary_2<Base, Arr_counting_traits_2<Base>>;
  using Counting_compare_x_near_boundary_2 =
    aos2::internal::Counting_compare_x_near_boundary_2<Base, Arr_counting_traits_2<Base>>;

  friend Counting_parameter_space_in_x_2;
  friend Counting_parameter_space_in_y_2;
  friend Counting_make_x_monotone_2;
  friend Counting_split_2;
  friend Counting_do_intersect_2;
  friend Counting_intersect_2;
  friend Counting_are_mergeable_2;
  friend Counting_merge_2;
  friend Counting_construct_opposite_2;
  friend Counting_construct_point_2;
  friend Counting_construct_x_monotone_curve_2;
  friend Counting_construct_curve_2;
  friend Counting_compare_endpoints_xy_2;
  friend Counting_approximate_2;
  friend Counting_is_on_x_identification_2;
  friend Counting_is_on_y_identification_2;
  friend Counting_compare_y_on_boundary_2;
  friend Counting_compare_y_near_boundary_2;
  friend Counting_x_on_boundary_2;
  friend Counting_compare_x_near_boundary_2;

public:
  enum Operation_id {
    COMPARE_X_2_OP = 0,
    COMPARE_XY_2_OP,
    CONSTRUCT_MIN_VERTEX_2_OP,
    CONSTRUCT_MAX_VERTEX_2_OP,
    IS_VERTICAL_2_OP,
    COMPARE_Y_AT_X_2_OP,
    EQUAL_2_POINTS_OP,
    EQUAL_2_CURVES_OP,
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
    APPROXIMATE_2_COORD_OP,
    APPROXIMATE_2_POINT_OP,
    APPROXIMATE_2_CURVE_OP,
    APPROXIMATE_2_BOUNDED_CURVE_OP,

    PARAMETER_SPACE_IN_X_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_X_2_POINT_OP,
    IS_ON_X_IDENTIFICATION_2_POINT_OP,
    IS_ON_X_IDENTIFICATION_2_CURVE_OP,
    COMPARE_Y_ON_BOUNDARY_2_OP,
    COMPARE_Y_NEAR_BOUNDARY_2_OP,

    PARAMETER_SPACE_IN_Y_2_CURVE_END_OP,
    PARAMETER_SPACE_IN_Y_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_POINT_OP,
    IS_ON_Y_IDENTIFICATION_2_CURVE_OP,
    COMPARE_X_ON_BOUNDARY_2_POINTS_OP,
    COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP,
    COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP,
    COMPARE_X_NEAR_BOUNDARY_2_OP,

    NUMBER_OF_OPERATIONS,
  };

  using Shared_base = std::shared_ptr<Base>;

  /*! constructs default.
   */
  template <typename ... Args>
  Arr_counting_traits_2(Args ... args) : m_base_traits(std::make_shared<Base>(std::forward<Args>(args)...)) {
    clear_counters();
    increment();
  }

  /*! constructs from a shared pointer.
   * \param[in] traits the taits being traced.
   * We use std::move to save an atomic operation.  Observe that moving a
   * shared_ptr simply transfers the internal pointer without touching the
   * atomic reference counter. This is faster than a copy.
   */
  Arr_counting_traits_2(Shared_base traits) : m_base_traits(std::move(traits)) {
    clear_counters();
    increment();
  }

  /*! disables copy constructor.
   */
  Arr_counting_traits_2(const Arr_counting_traits_2&) = delete;

  /*! obtains the counter of the given operation
   */
  std::size_t count(Operation_id id) const { return m_counters[id]; }

  /*! prints the counter associated with an operation.
   */
  template <typename OutStream>
  OutStream& print(OutStream& os, Operation_id id) const {
    if (! m_exist[id]) return os;
    os << m_names[id] << ": " << m_counters[id] << std::endl;
    return os;
  }

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

  /*! A functor that compares the \f$x\f$-coordinates of two points */
  class Compare_x_2 {
  private:
    typename Base::Compare_x_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_x_2(const Base& base, std::size_t& counter) : m_object(base.compare_x_2_object()), m_counter(counter) {}

    /*! operates
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      ++m_counter;
      return m_object(p1, p2);
    }
  };

  //! A functor that compares two points lexigoraphically: by \f$x\f$, then by \f$y\f$.
  class Compare_xy_2 {
  private:
    typename Base::Compare_xy_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_xy_2(const Base& base, std::size_t& counter) : m_object(base.compare_xy_2_object()), m_counter(counter) {}

    /*! operates
     */
    Comparison_result operator()(const Point_2& p1, const Point_2& p2) const {
      ++m_counter;
      return m_object(p1, p2);
    }
  };

  //! A functor that obtains the left endpoint of an \f$x\f$-monotone curve.
  class Construct_min_vertex_2 {
  private:
    typename Base::Construct_min_vertex_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Construct_min_vertex_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_min_vertex_2_object()), m_counter(counter) {}

    using Subcurve_ctr_minv = typename Base::Construct_min_vertex_2;

    /*! operates
     */
    decltype(std::declval<Subcurve_ctr_minv>().operator()(std::declval<X_monotone_curve_2>()))
    operator()(const X_monotone_curve_2& xcv) const {
      ++m_counter;
      return m_object(xcv);
    }
  };

  /*! A functor that obtains the right endpoint of an \f$x\f$-monotone curve. */
  class Construct_max_vertex_2 {
  private:
    typename Base::Construct_max_vertex_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Construct_max_vertex_2(const Base& base, std::size_t& counter) :
      m_object(base.construct_max_vertex_2_object()), m_counter(counter) {}

    using Subcurve_ctr_maxv = typename Base::Construct_max_vertex_2;

    /*! operates
     */
    decltype(std::declval<Subcurve_ctr_maxv>().operator()(std::declval<X_monotone_curve_2>()))
    operator()(const X_monotone_curve_2& xcv) const {
      ++m_counter;
      return m_object(xcv);
    }
  };

  /*! A functor that checks whether a given \f$x\f$-monotone curve is vertical.
   */
  class Is_vertical_2 {
  private:
    typename Base::Is_vertical_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Is_vertical_2(const Base& base, std::size_t& counter) : m_object(base.is_vertical_2_object()), m_counter(counter) {}

    /*! operates
     */
    bool operator()(const X_monotone_curve_2& xc) const {
      ++m_counter;
      return m_object(xc);
    }
  };

  /*! A functor that compares the \f$y\f$-coordinates of a point and an
   * \f$x\f$-monotone curve at the point \f$x\f$-coordinate.
   */
  class Compare_y_at_x_2 {
  private:
    typename Base::Compare_y_at_x_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_y_at_x_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_2_object()), m_counter(counter) {}

    /*! operates
     */
    Comparison_result operator()(const Point_2& p, const X_monotone_curve_2& xc) const {
      ++m_counter;
      return m_object(p, xc);
    }
  };

  /*! A functor that checks whether two points and two \f$x\f$-monotone curves
   * are identical.
   */
  class Equal_2 {
  private:
    typename Base::Equal_2 m_object;
    std::size_t& m_counter1;
    std::size_t& m_counter2;

  public:
    /*! constructs */
    Equal_2(const Base& base, std::size_t& counter1, std::size_t& counter2) :
      m_object(base.equal_2_object()), m_counter1(counter1), m_counter2(counter2) {}

    /*! operates
     */
    bool operator()(const Point_2& p1, const Point_2& p2) const {
      ++m_counter1;
      return m_object(p1, p2);
    }

    /*! operates
     */
    bool operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2) const {
      ++m_counter2;
      return m_object(xcv1, xcv2);
    }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the left of their intersection
   * point.
   */
  class Compare_y_at_x_left_2 {
  private:
    typename Base::Compare_y_at_x_left_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_y_at_x_left_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_left_2_object()), m_counter(counter) {}

    /*! operates
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const {
      ++m_counter;
      return m_object(xcv1, xcv2, p);
    }
  };

  /*! A functor that compares compares the \f$y\f$-coordinates of two
   * \f$x\f$-monotone curves immediately to the right of their intersection
   * point.
   */
  class Compare_y_at_x_right_2 {
  private:
    typename Base::Compare_y_at_x_right_2 m_object;
    std::size_t& m_counter;

  public:
    /*! constructs
     */
    Compare_y_at_x_right_2(const Base& base, std::size_t& counter) :
      m_object(base.compare_y_at_x_right_2_object()), m_counter(counter) {}

    /*! operates
     */
    Comparison_result operator()(const X_monotone_curve_2& xcv1, const X_monotone_curve_2& xcv2,
                                 const Point_2& p) const {
      ++m_counter;
      return m_object(xcv1, xcv2, p);
    }
  };

  Compare_x_2 compare_x_2_object() const
  { return Compare_x_2(traits(), m_counters[COMPARE_X_2_OP]); }

  Compare_xy_2 compare_xy_2_object() const
  { return Compare_xy_2(traits(), m_counters[COMPARE_XY_2_OP]); }

  Construct_min_vertex_2 construct_min_vertex_2_object() const
  { return Construct_min_vertex_2(traits(), m_counters[CONSTRUCT_MIN_VERTEX_2_OP]); }

  Construct_max_vertex_2 construct_max_vertex_2_object() const
  { return Construct_max_vertex_2(traits(), m_counters[CONSTRUCT_MAX_VERTEX_2_OP]); }

  Is_vertical_2 is_vertical_2_object() const
  { return Is_vertical_2(traits(), m_counters[IS_VERTICAL_2_OP]); }

  Compare_y_at_x_2 compare_y_at_x_2_object() const
  { return Compare_y_at_x_2(traits(), m_counters[COMPARE_Y_AT_X_2_OP]); }

  Equal_2 equal_2_object() const
  { return Equal_2(traits(), m_counters[EQUAL_2_POINTS_OP], m_counters[EQUAL_2_CURVES_OP]); }

  Compare_y_at_x_left_2 compare_y_at_x_left_2_object() const
  { return Compare_y_at_x_left_2(traits(), m_counters[COMPARE_Y_AT_X_LEFT_2_OP]); }

  Compare_y_at_x_right_2 compare_y_at_x_right_2_object() const
  { return Compare_y_at_x_right_2(traits(), m_counters[COMPARE_Y_AT_X_RIGHT_2_OP]); }

  /*! increments the construction counter.
   * \param doit indicates whether to actually increment the counter or not.
   * \return the counter at the end of the operation.
   */
  static std::size_t increment(bool doit = true) {
#ifdef CGAL_NO_ATOMIC
    static std::size_t counter;
#else
    static std::atomic<std::size_t> counter;
#endif
    if (doit) ++counter;
    return counter;
  }

  /*! cleans all operation counters.
   */
  void clear_counters() { m_counters = {}; }

  /*! obtains a const reference to the traits being counted.
   */
  const Base& traits() const { return *m_base_traits; }

  /*! obtains a reference to the traits being counted.
   */
  Base& traits() { return *m_base_traits; }

  /*! obtains the smart pointer to the traits being counted.
   */
  Shared_base shared_traits() const { return m_base_traits; }

private:
  //! The operation counters
  mutable std::array<std::size_t, NUMBER_OF_OPERATIONS> m_counters;
  const std::array<std::string, NUMBER_OF_OPERATIONS> m_names = {
    "COMPARE_X_2_OP",
    "COMPARE_XY_2_OP",
    "CONSTRUCT_MIN_VERTEX_2_OP",
    "CONSTRUCT_MAX_VERTEX_2_OP",
    "IS_VERTICAL_2_OP",
    "COMPARE_Y_AT_X_2_OP",
    "EQUAL_2_POINTS_OP",
    "EQUAL_2_CURVES_OP",
    "COMPARE_Y_AT_X_LEFT_2_OP",
    "COMPARE_Y_AT_X_RIGHT_2_OP",
    "MAKE_X_MONOTONE_2_OP",
    "SPLIT_2_OP",
    "DO_INTERSECT_2_OP",
    "INTERSECT_2_OP",
    "ARE_MERGEABLE_2_OP",
    "MERGE_2_OP",
    "CONSTRUCT_2_OPPOSITE_2_OP",
    "CONSTRUCT_POINT_2_OP",
    "CONSTRUCT_POINT_2_XY_OP",
    "CONSTRUCT_X_MONOTONE_CURVE_2_OP",
    "CONSTRUCT_CURVE_2_OP",
    "COMPARE_ENDPOINTS_XY_2_OP",
    "APPROXIMATE_2_COORD_OP",
    "APPROXIMATE_2_POINT_OP",
    "APPROXIMATE_2_CURVE_OP",
    "APPROXIMATE_2_BOUNDED_CURVE_OP",

    "PARAMETER_SPACE_IN_X_2_CURVE_END_OP",
    "PARAMETER_SPACE_IN_X_2_POINT_OP",
    "IS_ON_X_IDENTIFICATION_POINT_2_OP",
    "IS_ON_X_IDENTIFICATION_CURVE_2_OP",
    "COMPARE_Y_ON_BOUNDARY_2_OP",
    "COMPARE_Y_NEAR_BOUNDARY_2_OP",

    "PARAMETER_SPACE_IN_Y_2_CURVE_END_OP",
    "PARAMETER_SPACE_IN_Y_2_POINT_OP",
    "IS_ON_Y_IDENTIFICATION_2_POINT_OP",
    "IS_ON_Y_IDENTIFICATION_2_CURVE_OP",
    "COMPARE_X_ON_BOUNDARY_2_POINTS_OP",
    "COMPARE_X_ON_BOUNDARY_2_POINT_CURVE_END_OP",
    "COMPARE_X_ON_BOUNDARY_2_CURVE_ENDS_OP",
    "COMPARE_X_NEAR_BOUNDARY_2_OP"
  };

  const std::array<bool, NUMBER_OF_OPERATIONS> m_exist = {
    has_compare_x_2<Base>::value,
    has_compare_xy_2<Base>::value,
    has_construct_min_vertex_2<Base>::value,
    has_construct_max_vertex_2<Base>::value,
    has_is_vertical_2<Base>::value,
    has_compare_y_at_x_2<Base>::value,
    has_equal_2<Base>::value, // points
    has_equal_2<Base>::value, // curves
    has_compare_y_at_x_left_2<Base>::value,
    has_compare_y_at_x_right_2<Base>::value,
    has_make_x_monotone_2<Base>::value,
    has_split_2<Base>::value,
    has_do_intersect_2<Base>::value,
    has_intersect_2<Base>::value,
    has_are_mergeable_2<Base>::value,
    has_merge_2<Base>::value,
    has_construct_opposite_2<Base>::value,
    has_construct_point_2<Base>::value,
    has_construct_point_2_xy<Base>::value,
    has_construct_x_monotone_curve_2<Base>::value,
    has_construct_curve_2<Base>::value,
    has_compare_endpoints_xy_2<Base>::value,
    has_approximate_2<Base>::value, // coordinate
    has_approximate_2_point<Base>::value, // point
    has_approximate_2_xcv<Base>::value, // curve
    has_approximate_2_xcv_bounds<Base>::value, // bounded curve
    has_parameter_space_in_x_2<Base>::value, // curve end
    has_parameter_space_in_x_2<Base>::value, // point
    has_is_on_x_identification_2<Base>::value, // point
    has_is_on_x_identification_2<Base>::value, // curve
    has_compare_y_on_boundary_2<Base>::value,
    has_compare_y_near_boundary_2<Base>::value,
    has_parameter_space_in_y_2<Base>::value, // curve end
    has_parameter_space_in_y_2<Base>::value, // point
    has_is_on_y_identification_2<Base>::value, // point
    has_is_on_y_identification_2<Base>::value, // curve
    has_compare_x_on_boundary_2<Base>::value, // points
    has_compare_x_on_boundary_2<Base>::value, // point, curve end
    has_compare_x_on_boundary_2<Base>::value, // curve ends
    has_compare_x_near_boundary_2<Base>::value
  };

  //! The traitse being counted
  Shared_base m_base_traits;
};

template <typename OutStream, class BaseTraits>
inline OutStream& operator<<(OutStream& os, const Arr_counting_traits_2<BaseTraits>& traits) {
  using Traits = Arr_counting_traits_2<BaseTraits>;
  std::size_t sum = 0;
  for (auto i = 0; i < Traits::NUMBER_OF_OPERATIONS; ++i) {
    sum += traits.count(static_cast<typename Traits::Operation_id>(i));
    traits.print(os, static_cast<typename Traits::Operation_id>(i));
  }
  os << "total # = " << sum << std::endl
     << "# of traits constructed = " << Traits::increment(false)
     << std::endl;
  return os;
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
