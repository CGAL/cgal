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
#ifndef GENERALIZED_MAP_FWD_H
#define GENERALIZED_MAP_FWD_H 1

#include <CGAL/memory.h>
#include <CGAL/tags.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

template<unsigned int d_, class Items_, class Alloc_>
class Generalized_map_storage_1;

template<unsigned int d_, class Items_, class Alloc_>
class Generalized_map_storage_with_index;

struct Generic_map_min_items;

namespace internal
{
template<typename Tag>
struct Default_storage_for_gmap_when_tag
{
  template<unsigned int d_, class Items_, class Alloc_>
  using type=Generalized_map_storage_1<d_, Items_, Alloc_>;
};
template<>
struct Default_storage_for_gmap_when_tag<CGAL::Tag_true>
{
  template<unsigned int d_, class Items_, class Alloc_>
  using type=Generalized_map_storage_with_index<d_, Items_, Alloc_>;
};

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_use_index_tag_gmap,Use_index,false)
template<typename T, bool typedefined=Has_use_index_tag_gmap<T>::value>
struct Default_storage_for_gmap
{
  template<unsigned int d_, class Items_, class Alloc_>
  using type=Generalized_map_storage_1<d_, Items_, Alloc_>;
};
template<typename T>
struct Default_storage_for_gmap<T, true>
{
  template<unsigned int d_, class Items_, class Alloc_>
  using type=typename CGAL::internal::template
  Default_storage_for_gmap_when_tag<typename T::Use_index>::
  template type<d_, Items_, Alloc_>;
};
} // namespace internal

template<unsigned int d_, class Refs_,
         class Items_=Generic_map_min_items,
         class Alloc_=CGAL_ALLOCATOR(int),
         class Storage_=typename internal::template
         Default_storage_for_gmap<Items_>::template type<d_, Items_, Alloc_>>
class Generalized_map_base;

template <unsigned int d_,
          class Items_=Generic_map_min_items,
          class Alloc_=CGAL_ALLOCATOR(int),
         class Storage_=typename internal::template
          Default_storage_for_gmap<Items_>::template type<d_, Items_, Alloc_>>
class Generalized_map;

} // CGAL

#endif // GENERALIZED_MAP_FWD_H
