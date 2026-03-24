// Copyright (c) 2026 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// author(s)     : Léo Valque

#ifndef CGAL_SNAP_ROUNDING_HELPERS_2_H
#define CGAL_SNAP_ROUNDING_HELPERS_2_H

#include <CGAL/license/Snap_rounding_2.h>
#include <CGAL/Polygon_2.h>

namespace CGAL{
namespace internal{

template <class T>
inline constexpr bool is_instance_of_Polygon_2 = false;

template <class K, class C>
inline constexpr bool is_instance_of_Polygon_2< CGAL::Polygon_2<K, C> > = true;

}
}

#endif