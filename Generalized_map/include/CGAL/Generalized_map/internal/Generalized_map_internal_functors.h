// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_GENERALIZED_MAP_INTERNAL_FUNCTORS_H
#define CGAL_GENERALIZED_MAP_INTERNAL_FUNCTORS_H

/* Definition of functors used internally for generalized maps.
 *
 * internal::Alpha_functor<Dart, i...> to call several beta on the given dart.
 *   Indices are given as parameter of the run function.
 *
 * internal::Alpha_functor_static<Dart, i...> to call several beta on the given
 *   dart. Indices are given as template arguments.
 *
 */

#include <CGAL/Combinatorial_map/internal/Combinatorial_map_internal_functors.h>

namespace CGAL
{
// ****************************************************************************
namespace internal
{
// ****************************************************************************
// Alpha functor, used to combine several alpha.
template<typename GMap, typename Dart_descriptor, typename ... Alphas>
struct Alpha_functor;

template<typename GMap, typename Dart_descriptor>
struct Alpha_functor<GMap, Dart_descriptor, int>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart, int B)
  { return AMap.get_alpha(ADart, B); }
};

template<typename GMap, typename Dart_descriptor>
struct Alpha_functor<GMap, Dart_descriptor, unsigned int>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart, unsigned int B)
  { return  AMap.get_alpha(ADart, B); }
};

template<typename GMap, typename Dart_descriptor, typename ... Alphas>
struct Alpha_functor<GMap, Dart_descriptor, int, Alphas...>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart, int B, Alphas... alphas)
  { return Alpha_functor<GMap, Dart_descriptor, Alphas...>::
      run(AMap, AMap.get_alpha(ADart, B), alphas...); }
};

template<typename GMap, typename Dart_descriptor, typename ... Alphas>
struct Alpha_functor<GMap, Dart_descriptor, unsigned int, Alphas...>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart, unsigned int B,
                         Alphas... alphas)
  { return Alpha_functor<GMap, Dart_descriptor, Alphas...>::
      run(AMap, AMap.get_alpha(ADart, B), alphas...); }
};
// ****************************************************************************
template<typename GMap, typename Dart_descriptor, int ... Alphas>
struct Alpha_functor_static;

template<typename GMap, typename Dart_descriptor, int B>
struct Alpha_functor_static<GMap, Dart_descriptor, B>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart)
  { return AMap.template get_alpha<B>(ADart); }
};

template<typename GMap, typename Dart_descriptor, int B, int ... Alphas>
struct Alpha_functor_static<GMap, Dart_descriptor, B, Alphas...>
{
  static Dart_descriptor run(GMap& AMap, Dart_descriptor ADart)
  { return Alpha_functor_static<GMap, Dart_descriptor, Alphas...>::
        run(AMap, AMap.template get_alpha<B>(ADart)); }
};
// ****************************************************************************
} // namespace internal
} // namespace CGAL

#endif // CGAL_GENERALIZED_MAP_INTERNAL_FUNCTORS_H
