// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_OCTREE_H

#include <CGAL/license/Shape_detection.h>
#include <CGAL/Octree.h>

namespace CGAL {
namespace Shape_detection {

template<class Traits>
class Octree : CGAL::Octree::Octree<typename Traits::Input_range, typename Traits::Point_map> {

};

}
}

#endif
