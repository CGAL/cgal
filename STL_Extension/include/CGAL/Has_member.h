// Copyright (c) 2017 Inria (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Clement Jamin

#ifndef CGAL_HAS_MEMBER_H
#define CGAL_HAS_MEMBER_H

// Macro used to check if a type T has a member named `X`
// It generates a class has_X<T> where has_X<T>::value is a boolean
// See example in Concurrent_compact_container.h
#define CGAL_GENERATE_MEMBER_DETECTOR(X)                                                           \
template<typename T, typename... Args>                                                             \
class has_##X                                                                                      \
{                                                                                                  \
  struct Fallback { int X; };                                                                      \
  struct Derived : T, Fallback { };                                                                \
                                                                                                   \
  template<typename U, U>                                                                          \
  struct Check;                                                                                    \
                                                                                                   \
  typedef char ArrayOfOne[1];                                                                      \
  typedef char ArrayOfTwo[2];                                                                      \
                                                                                                   \
  template<typename U>                                                                             \
  static ArrayOfOne& func(Check<int Fallback::*, &U::X> *);                                        \
                                                                                                   \
  template<typename U> static ArrayOfTwo & func(...);                                              \
                                                                                                   \
public:                                                                                            \
    typedef has_##X type;                                                                          \
    enum { value = sizeof(func<Derived>(0)) == 2 };                                                \
};                                                                                                 \
                                                                                                   \
template<typename C, typename Ret, typename... Args>                                               \
struct has_##X<C, Ret(Args...)>                                                                    \
{                                                                                                  \
private:                                                                                           \
  template<typename T>                                                                             \
  static constexpr auto check(T*)                                                                  \
    -> typename std::is_same<decltype(std::declval<T>().X(std::declval<Args>()...)), Ret>::type;   \
                                                                                                   \
  template<typename>                                                                               \
  static constexpr std::false_type check(...);                                                     \
                                                                                                   \
  typedef decltype(check<C>(0)) type;                                                              \
                                                                                                   \
public:                                                                                            \
  static constexpr bool value = type::value;                                                       \
};

#endif // CGAL_HAS_MEMBER_H
