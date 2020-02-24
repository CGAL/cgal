// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_2_H
#define CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_2_H

#include <CGAL/license/Alpha_shapes_2.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/Weighted_alpha_shape_euclidean_traits_2.h>"
#define CGAL_DEPRECATED_MESSAGE_DETAILS \
  "The kernel K can be used directly as traits since weighted points and "\
  "the associated function objects are now part of the concept Kernel."
#include <CGAL/internal/deprecation_warning.h>

namespace CGAL {

template< class K_ >
class Weighted_alpha_shape_euclidean_traits_2
  : public K_
{
public:
  Weighted_alpha_shape_euclidean_traits_2() { }
  Weighted_alpha_shape_euclidean_traits_2(const K_& k) : K_(k) { }
};

} // namespace CGAL

#endif // CGAL_WEIGHTED_ALPHA_SHAPE_EUCLIDEAN_TRAITS_2_H
