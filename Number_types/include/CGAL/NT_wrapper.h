// Copyright (c) 2024, GeometryFactory Sarl (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_NT_WRAPPER_H
#define CGAL_NT_WRAPPER_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/functional.h>
#include <functional>
#include <source_location>

namespace CGAL {

/*
 * This class template `NT_wapper<NT>` is currently undocumented, on purpose.
 *
 * This class template provides a wrapper for number types. It calls a function
 * `f` when a constructor, destructor, or assignment operator is called. This
 * is useful to detect when a copy constructor is called, for example. An
 * example of use case is in the test file
 * `Number_types/test/Number_types/wrapping_type.cpp`.
 */

template <typename NT>
class NT_wrapper {
  NT value;

  static void call_f(const std::source_location& l = std::source_location::current()) {
    if (f)
      f(l);
  }

public:
  inline static std::function<void(const std::source_location&)> f;

  NT_wrapper() : value(0) { call_f(std::source_location::current()); }
  NT_wrapper(const NT& val) : value(val) { call_f(std::source_location::current());}
  NT_wrapper(NT&& val) : value(std::move(val)) { call_f(std::source_location::current()); }
  NT_wrapper(int val) : value(val) { call_f(std::source_location::current()); }

  template <typename T, typename = std::enable_if_t<!std::is_same_v<T, NT>>>
  NT_wrapper(const NT_wrapper<T>& other) : value(other.get_value()) { call_f(std::source_location::current()); }

  NT_wrapper(const NT_wrapper& other) : value(other.value) { call_f(std::source_location::current()); }
  NT_wrapper(NT_wrapper&& other) : value(std::move(other.value)) { call_f(std::source_location::current()); }
  ~NT_wrapper() { call_f(std::source_location::current()); }

  NT_wrapper& operator=(const NT_wrapper& other) {
    if (this != &other) {
      value = other.value;
      call_f(std::source_location::current());
    }
    return *this;
  }

  NT_wrapper& operator=(NT_wrapper&& other) {
    if (this != &other) {
      value = std::move(other.value);
      call_f(std::source_location::current());
    }
    return *this;
  }

  NT get_value() const {
    return value;
  }

  operator NT() const {
    return value;
  }

  NT_wrapper operator+(const NT_wrapper& rhs) const {
    return NT_wrapper(value + rhs.value);
  }

  NT_wrapper operator-(const NT_wrapper& rhs) const {
    return NT_wrapper(value - rhs.value);
  }

  NT_wrapper operator*(const NT_wrapper& rhs) const {
    return NT_wrapper(value * rhs.value);
  }

  NT_wrapper operator/(const NT_wrapper& rhs) const {
    return NT_wrapper(value / rhs.value);
  }

  NT_wrapper operator/(int rhs) const {
    return NT_wrapper(value / rhs);
  }

  bool operator==(const NT_wrapper& rhs) const {
    return value == rhs.value;
  }

  bool operator!=(const NT_wrapper& rhs) const {
    return value != rhs.value;
  }

  bool operator<(const NT_wrapper& rhs) const {
    return value < rhs.value;
  }

  bool operator>(const NT_wrapper& rhs) const {
    return value > rhs.value;
  }

  bool operator<=(const NT_wrapper& rhs) const {
    return value <= rhs.value;
  }

  bool operator>=(const NT_wrapper& rhs) const {
    return value >= rhs.value;
  }
};

CGAL_DEFINE_COERCION_TRAITS_FROM_TO_TEM(NT_wrapper<NT>, double, typename NT)

template <typename NT>
struct Real_embeddable_traits<NT_wrapper<NT>>
    : public INTERN_RET::Real_embeddable_traits_base<NT_wrapper<NT>,
                                                     typename Real_embeddable_traits<NT>::Is_real_embeddable>
{
  using Type = NT_wrapper<NT>;
  class Sgn
    : public CGAL::cpp98::unary_function< NT, ::CGAL::Sign > {
  public:
    ::CGAL::Sign operator()( const Type& x ) const {
      return CGAL_NTS sign( x.get_value() );
    }
  };

  class To_double
    : public CGAL::cpp98::unary_function< NT, double > {
  public:
    double operator()( const Type& x) const {
      return CGAL_NTS to_double( x.get_value() );
    }
  };

  class To_interval
    : public CGAL::cpp98::unary_function< NT, std::pair<double,double> > {
  public:
    std::pair<double,double> operator()( const Type& x ) const {
      return CGAL_NTS to_interval( x.get_value() );
    }
  };
};

namespace internal {

template <typename T>
struct Exact_field_selector<NT_wrapper<T>> {
  using Type = NT_wrapper<typename Exact_field_selector<T>::Type>;
  using type = Type;
};

template <typename T>
struct Exact_ring_selector<NT_wrapper<T>> {
  using Type = NT_wrapper<typename Exact_ring_selector<T>::Type>;
  using type = Type;
};

} // namespace internal

struct Wrapped_epick
    : public CGAL::internal::Epick_with_filtered_predicates<NT_wrapper<double>, Wrapped_epick>
{};
struct Wrapped_epeck
    : public CGAL::internal::Epeck_lazy_base_with_type_equality<NT_wrapper<CGAL::Epeck_ft>, Wrapped_epeck>
{};

} // namespace CGAL

#endif // CGAL_NT_WRAPPER_H