// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_UNCERTAIN_H
#define CGAL_UNCERTAIN_H

#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <stdexcept>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  // Accessory traits class to provide min and max value of a type.
  // Specialized for bool and CGAL's enums.
  template < typename T >
  struct Minmax_traits;

  template <>
  struct Minmax_traits <bool>
  {
    static const bool min = false;
    static const bool max = true;
  };

  template <>
  struct Minmax_traits <Sign>
  {
    static const Sign min = NEGATIVE;
    static const Sign max = POSITIVE;
  };

  template <>
  struct Minmax_traits <Bounded_side>
  {
    static const Bounded_side min = ON_UNBOUNDED_SIDE;
    static const Bounded_side max = ON_BOUNDED_SIDE;
  };

  template <>
  struct Minmax_traits <Angle>
  {
    static const Angle min = OBTUSE;
    static const Angle max = ACUTE;
  };

} // namespace CGALi


// Exception type for the automatic conversion.
class Uncertain_conversion_exception
  : std::range_error
{
public:
  Uncertain_conversion_exception(const std::string &s)
    : std::range_error(s) {}

  ~Uncertain_conversion_exception() throw() {}
};


// Encodes an interval [min;max] of values of type T.
// The primary template is supposed to work for enums and bool.

template < typename T >
class Uncertain
{
  T _i, _s;

  static unsigned failures;   // Number of conversion failures.

public:

  typedef CGAL::Uncertain_conversion_exception  Uncertain_conversion_exception;

  typedef T  value_type;

  Uncertain()
    : _i(CGALi::Minmax_traits<T>::min),
      _s(CGALi::Minmax_traits<T>::max) {}

  Uncertain(T t)
    : _i(t), _s(t) {}

  Uncertain(T i, T s)
    : _i(i), _s(s) {}


  T inf() const { return _i; }
  T sup() const { return _s; }

  bool is_certain() const { return _i == _s; }

  T make_certain() const
  {
    if (is_certain())
      return _i;
    ++Uncertain<bool>::number_of_failures();
    throw Uncertain_conversion_exception(
                        "undecidable conversion of CGAL::Uncertain<T>");
  }

#if 1 // Comment out in order to spot some unwanted conversions.
  // NB : the general conversion to bool might be too risky.
  //      boost::tribool uses a more restricted conversion (see below).
  operator T() const
  {
    return make_certain();
  }
#else
private:
  struct dummy {
    void nonnull() {};
  };

  typedef void (dummy::*safe_bool)();

public:
  operator safe_bool() const
  {
    return make_certain() ? &dummy::nonnull : 0;
  }
#endif

  static unsigned & number_of_failures() { return failures; }

  static Uncertain indeterminate();
};


template < typename T >
unsigned
Uncertain<T>::failures = 0;


// Access functions
// ----------------

template < typename T >
inline
T inf(Uncertain<T> i)
{
  return i.inf();
}

template < typename T >
inline
T sup(Uncertain<T> i)
{
  return i.sup();
}


// possibly() declarations (needed as used in the assertions)
// ----------------------------------------------------------

inline bool possibly(bool b);
inline bool possibly(Uncertain<bool> c);


// Basic functions
// ---------------

#if 1
// Commenting them out allows the test-suite to spot those predicates that do not
// propagate Uncertain-ty correctly.  But they should probably be enabled
// for now, for backward-compatibility.
template < typename T >
inline
bool is_certain(T)
{
  return true;
}

template < typename T >
inline
T get_certain(T t)
{
  return t;
}
#endif

template < typename T >
inline
bool is_certain(Uncertain<T> a)
{
  return a.is_certain();
}

template < typename T >
inline
T get_certain(Uncertain<T> a)
{
  CGAL_assertion(is_certain(a));
  return a.inf();
}

template < typename T >
inline
Uncertain<T>
Uncertain<T>::indeterminate()
{
  return Uncertain<T>();
}

template < typename T >
inline
bool is_indeterminate(T)
{
  return false;
}

template < typename T >
inline
bool is_indeterminate(Uncertain<T> a)
{
  return ! a.is_certain();
}


// certainly/possibly
// ------------------

inline bool certainly(bool b) { return b; }
inline bool possibly(bool b) { return b; }

inline
bool certainly(Uncertain<bool> c)
{
  if (is_certain(c))
    return get_certain(c);
  return false;
}

inline
bool possibly(Uncertain<bool> c)
{
  if (is_certain(c))
    return get_certain(c);
  return true;
}


// Boolean operations for Uncertain<bool>
// --------------------------------------

inline
Uncertain<bool> operator!(Uncertain<bool> a)
{
  return Uncertain<bool>(!a.sup(), !a.inf());
}

inline
Uncertain<bool> operator|(Uncertain<bool> a, Uncertain<bool> b)
{
  return Uncertain<bool>(a.inf() || b.inf(), a.sup() || b.sup());
}

inline
Uncertain<bool> operator|(bool a, Uncertain<bool> b)
{
  return Uncertain<bool>(a || b.inf(), a || b.sup());
}

inline
Uncertain<bool> operator|(Uncertain<bool> a, bool b)
{
  return Uncertain<bool>(a.inf() || b, a.sup() || b);
}

inline
Uncertain<bool> operator&(Uncertain<bool> a, Uncertain<bool> b)
{
  return Uncertain<bool>(a.inf() && b.inf(), a.sup() && b.sup());
}

inline
Uncertain<bool> operator&(bool a, Uncertain<bool> b)
{
  return Uncertain<bool>(a && b.inf(), a && b.sup());
}

inline
Uncertain<bool> operator&(Uncertain<bool> a, bool b)
{
  return Uncertain<bool>(a.inf() && b, a.sup() && b);
}


// Equality operators

template < typename T >
inline
Uncertain<bool> operator==(Uncertain<T> a, Uncertain<T> b)
{
  if (is_indeterminate(a) || is_indeterminate(b))
    return Uncertain<bool>::indeterminate();
  return a.inf() == b.inf();
}

template < typename T >
inline
Uncertain<bool> operator==(Uncertain<T> a, T b)
{
  if (is_indeterminate(a))
    return Uncertain<bool>::indeterminate();
  return a.inf() == b;
}

template < typename T >
inline
Uncertain<bool> operator==(T a, Uncertain<T> b)
{
  if (is_indeterminate(b))
    return Uncertain<bool>::indeterminate();
  return a == b.inf();
}

template < typename T >
inline
Uncertain<bool> operator!=(Uncertain<T> a, Uncertain<T> b)
{
  return ! (a == b);
}

template < typename T >
inline
Uncertain<bool> operator!=(Uncertain<T> a, T b)
{
  return ! (a == b);
}

template < typename T >
inline
Uncertain<bool> operator!=(T a, Uncertain<T> b)
{
  return ! (a == b);
}


// Comparison operators (useful for enums only, I guess (?) ).

template < typename T >
inline
Uncertain<bool> operator<(Uncertain<T> a, Uncertain<T> b)
{
  if (a.sup() < b.inf())
    return true;
  if (a.inf() >= b.sup())
    return false;
  return Uncertain<bool>::indeterminate();
}

template < typename T >
inline
Uncertain<bool> operator<(Uncertain<T> a, T b)
{
  if (a.sup() < b)
    return true;
  if (a.inf() >= b)
    return false;
  return Uncertain<bool>::indeterminate();
}

template < typename T >
inline
Uncertain<bool> operator<(T a, Uncertain<T> b)
{
  if (a < b.inf())
    return true;
  if (a >= b.sup())
    return false;
  return Uncertain<bool>::indeterminate();
}

template < typename T >
inline
Uncertain<bool> operator>(Uncertain<T> a, Uncertain<T> b)
{
  return b < a;
}

template < typename T >
inline
Uncertain<bool> operator>(Uncertain<T> a, T b)
{
  return b < a;
}

template < typename T >
inline
Uncertain<bool> operator>(T a, Uncertain<T> b)
{
  return b < a;
}

template < typename T >
inline
Uncertain<bool> operator<=(Uncertain<T> a, Uncertain<T> b)
{
  return !(b < a);
}

template < typename T >
inline
Uncertain<bool> operator<=(Uncertain<T> a, T b)
{
  return !(b < a);
}

template < typename T >
inline
Uncertain<bool> operator<=(T a, Uncertain<T> b)
{
  return !(b < a);
}

template < typename T >
inline
Uncertain<bool> operator>=(Uncertain<T> a, Uncertain<T> b)
{
  return !(a < b);
}

template < typename T >
inline
Uncertain<bool> operator>=(Uncertain<T> a, T b)
{
  return !(a < b);
}

template < typename T >
inline
Uncertain<bool> operator>=(T a, Uncertain<T> b)
{
  return !(a < b);
}


// Maker function (a la std::make_pair).

template < typename T >
inline
Uncertain<T> make_uncertain(T t)
{
  return Uncertain<T>(t);
}

template < typename T >
inline
Uncertain<T> make_uncertain(Uncertain<T> t)
{
  return t;
}


// make_certain() : Forcing a cast to certain (possibly throwing).
// This is meant to be used only in cases where we cannot easily propagate the
// uncertainty, such as when used in a switch statement (code may later be
// revisited to do things in a better way).

template < typename T >
inline
T make_certain(T t)
{
  return t;
}

template < typename T >
inline
T make_certain(Uncertain<T> t)
{
  return t.make_certain();
}


// opposite
template < typename T > // should be constrained only for enums.
inline
Uncertain<T> operator-(Uncertain<T> u)
{
  return Uncertain<T>(-u.sup(), -u.inf());
}

// "sign" multiplication.
// Should be constrained only for "sign" enums, useless for bool.
template < typename T >
Uncertain<T> operator*(Uncertain<T> a, Uncertain<T> b)
{
  if (a.inf() >= 0)                                   // e>=0
  {
    // b>=0     [a.inf()*b.inf(); a.sup()*b.sup()]
    // b<=0     [a.sup()*b.inf(); a.inf()*b.sup()]
    // b~=0     [a.sup()*b.inf(); a.sup()*b.sup()]
    T aa = a.inf(), bb = a.sup();
    if (b.inf() < 0)
    {
        aa = bb;
        if (b.sup() < 0)
            bb = a.inf();
    }
    return Uncertain<T>(aa * b.inf(), bb * b.sup());
  }
  else if (a.sup()<=0)                                // e<=0
  {
    // b>=0     [a.inf()*b.sup(); a.sup()*b.inf()]
    // b<=0     [a.sup()*b.sup(); a.inf()*b.inf()]
    // b~=0     [a.inf()*b.sup(); a.inf()*b.inf()]
    T aa = a.sup(), bb = a.inf();
    if (b.inf() < 0)
    {
        aa=bb;
        if (b.sup() < 0)
            bb=a.sup();
    }
    return Uncertain<T>(bb * b.sup(), aa * b.inf());
  }
  else                                          // 0 \in [inf();sup()]
  {
    if (b.inf()>=0)                           // d>=0
      return Uncertain<T>(a.inf() * b.sup(), a.sup() * b.sup());
    if (b.sup()<=0)                           // d<=0
      return Uncertain<T>(a.sup() * b.inf(), a.inf() * b.inf());
                                                // 0 \in d
    T tmp1 = a.inf() * b.sup();
    T tmp2 = a.sup() * b.inf();
    T tmp3 = a.inf() * b.inf();
    T tmp4 = a.sup() * b.sup();
    return Uncertain<T>((std::min)(tmp1, tmp2), (std::max)(tmp3, tmp4));
  }
}

template < typename T >
inline
Uncertain<T> operator*(T a, Uncertain<T> b)
{
	return Uncertain<T>(a) * b;
}

template < typename T >
inline
Uncertain<T> operator*(Uncertain<T> a, T b)
{
	return a * Uncertain<T>(b);
}

// enum_cast overload

#ifdef CGAL_CFG_MATCHING_BUG_5

template < typename T, typename U >
inline
Uncertain<T> enum_cast_bug(Uncertain<U> u, const T*)
{
  return Uncertain<T>(static_cast<const T>(u.inf()),
                      static_cast<const T>(u.sup()));
}

#else

template < typename T, typename U >
inline
Uncertain<T> enum_cast(Uncertain<U> u)
{
  return Uncertain<T>(static_cast<T>(u.inf()), static_cast<T>(u.sup()));
}

#endif

CGAL_END_NAMESPACE

#endif // CGAL_UNCERTAIN_H
