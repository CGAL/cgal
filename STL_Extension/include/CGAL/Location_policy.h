// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_LOCATION_POLICY_H
#define CGAL_LOCATION_POLICY_H

#include <CGAL/Complexity_tags.h>

namespace CGAL {

// A policy to select the complexity of point location of a data-structure.

template < typename Tag >
struct Location_policy {};

typedef Location_policy<Fast>     Fast_location;
typedef Location_policy<Compact>  Compact_location;

} // namespace CGAL

#endif // CGAL_LOCATION_POLICY_H
