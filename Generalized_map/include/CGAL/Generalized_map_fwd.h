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

namespace CGAL {

template<unsigned int d_, class Items_, class Alloc_ >
class Generalized_map_storage_1;

struct Generic_map_min_items;

template < unsigned int d_, class Refs,
           class Items_=Generic_map_min_items,
           class Alloc_=CGAL_ALLOCATOR(int),
           class Storage_= Generalized_map_storage_1<d_, Items_, Alloc_> >
class Generalized_map_base;

template < unsigned int d_,
           class Items_=Generic_map_min_items,
           class Alloc_=CGAL_ALLOCATOR(int),
           class Storage_= Generalized_map_storage_1<d_, Items_, Alloc_> >
class Generalized_map;

} // CGAL

#endif // GENERALIZED_MAP_FWD_H
