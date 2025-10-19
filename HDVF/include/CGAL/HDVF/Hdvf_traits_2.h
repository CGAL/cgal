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

#ifndef CGAL_HDVF_TRAITS_2_H
#define CGAL_HDVF_TRAITS_2_H

#include <CGAL/license/HDVF.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Dimension.h>
#include <vector>

namespace CGAL {
namespace Homological_discrete_vector_field {

  /*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_traits_2` implements the `HDVFTraits` concept for 2D data, using a geometric kernel `K`.

 @tparam K a geometric kernel model of the `Kernel` concept.

 \cgalModels{HDVFTraits}

 */

template <typename K>
struct Hdvf_traits_2 {
    using Dimension = Dimension_tag< 2 >;
    typedef K Kernel;
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;
    typedef typename K::FT FT;
    typedef CGAL::Bbox_2 Bbox;

    Point to_point(const std::vector<FT>& coords) const {
        return Point(coords[0], coords[1]);
    }
};

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_TRAITS_2_H
