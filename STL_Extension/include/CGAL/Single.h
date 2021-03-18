// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_SINGLE_H
#define CGAL_SINGLE_H

namespace CGAL
{

// Class to create a range of one single element
template <typename T>
class Single
{
public:

  typedef T* iterator;
  typedef T* const_iterator;
  typedef std::size_t size_type;

private:

  T& t;

public:

  Single (T& t) : t(t) { }

  iterator begin() const { return &t; }
  iterator end()   const { return &t + 1; }

  constexpr size_type size()  const { return 1; }
  constexpr bool      empty() const { return false; }
};

template <typename T>
Single<const T> make_single (const T& t)
{
  return Single<const T>(t);
}

template <typename T>
Single<T> make_single (T& t)
{
  return Single<T>(t);
}



}

#endif // CGAL_SINGLE_H
