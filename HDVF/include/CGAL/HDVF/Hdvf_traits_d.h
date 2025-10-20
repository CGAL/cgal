// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_TRAITS_D_H
#define CGAL_HDVF_TRAITS_D_H

#include <CGAL/license/HDVF.h>
#include <CGAL/Bbox_d.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Homological_discrete_vector_field {

  /*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_traits_d` implements the `HDVFTraits` concept for dD data, using a geometric kernel `K`.

 @tparam K a geometric traits class. Must be either `Epick_d` or `Epeck_d` with a fixed dimension.

 \cgalModels{HDVFTraits}

 */

template <typename K>
struct Hdvf_traits_d {
    using Dimension = Dimension_tag< K::Dimension::value >;
    typedef K Kernel;
    typedef typename K::Point_d Point;
    typedef typename K::FT FT;
    typedef CGAL::Bbox_d<Dimension> Bbox;
    typedef typename K::Point_3 Point3;
    static Point3 to_point3(const Point& p) { return Point3(p[0], p[1], p[2]); }
    }
};

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_TRAITS_D_H
