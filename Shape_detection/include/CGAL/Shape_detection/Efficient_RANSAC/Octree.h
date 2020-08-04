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

#include <stack>
#include <limits>

#include <CGAL/Random.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>

#include <CGAL/Octree.h>
#include <CGAL/Octree/Node.h>

namespace CGAL {
namespace Shape_detection {
namespace internal {


template<class Traits>
class Direct_octree : public Octree::Octree<typename Traits::Input_range, typename Traits::Point_map> {

  std::size_t m_offset;

public:

  std::size_t size() const {
    return this->root().size();
  }

  const Bbox_3 &boundingBox() const {
    return this->bbox(this->root());
  }

  std::size_t maxLevel() const {
    return this->max_depth_reached();
  }

  std::size_t offset() const { return m_offset; }
};

template<class Traits>
class Indexed_octree : public Octree::Octree<typename Traits::Input_range, typename Traits::Point_map> {

public:

  std::size_t size() const {
    return this->root().size();
  }

  const Bbox_3 &boundingBox() const {
    return this->bbox(this->root());
  }

  std::size_t maxLevel() const {
    return this->max_depth_reached();
  }
};

}
}
}

#endif
