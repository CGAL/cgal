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

#ifndef CGAL_HDVF_TRAITS_3_H
#define CGAL_HDVF_TRAITS_3_H

#include <functional>
#include <CGAL/license/HDVF.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Dimension.h>
#include <vector>

namespace CGAL {
namespace Homological_discrete_vector_field {

  /*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Hdvf_traits_3` implements the `HDVFTraits` concept for 3D data, using a geometric kernel `K`.

 @tparam K a geometric kernel model of the `Kernel` concept.

 \cgalModels{HDVFTraits}

 */

template <typename K>
struct Hdvf_traits_3 {
    using Dimension = Dimension_tag< 3 >;
    typedef K Kernel;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::FT FT;
    typedef CGAL::Bbox_3 Bbox;
    typedef typename K::Point_3 Point3;
    static std::function<Point3(const Point&)> to_point3;
};

template <typename K>
std::function<typename K::Point_3(const typename K::Point_3&)> Hdvf_traits_3<K>::to_point3 = [](const K::Point_3& p) { return typename K::Point_3(p); };

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */

#endif // CGAL_HDVF_TRAITS_3_H
