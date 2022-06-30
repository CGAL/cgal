// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_FEATURES_DETECTION_CUBE_ISOMETRIES_H
#define CGAL_MESH_3_FEATURES_DETECTION_CUBE_ISOMETRIES_H

#include <CGAL/license/Mesh_3.h>

#include <array>

namespace CGAL
{
namespace Mesh_3
{
namespace internal
{
  using Permutation = std::array<std::uint8_t, 8>;

  constexpr Permutation cube_isometries[] = {
      {0,1,2,3,4,5,6,7},
      {1,0,3,2,5,4,7,6},
      {4,5,0,1,6,7,2,3},
      {5,4,1,0,7,6,3,2},
      {6,7,4,5,2,3,0,1},
      {7,6,5,4,3,2,1,0},
      {2,3,6,7,0,1,4,5},
      {3,2,7,6,1,0,5,4},
      {1,5,3,7,0,4,2,6},
      {5,1,7,3,4,0,6,2},
      {5,4,7,6,1,0,3,2},
      {4,5,6,7,0,1,2,3},
      {4,0,6,2,5,1,7,3},
      {0,4,2,6,1,5,3,7},
      {1,3,0,2,5,7,4,6},
      {3,1,2,0,7,5,6,4},
      {3,2,1,0,7,6,5,4},
      {2,3,0,1,6,7,4,5},
      {2,0,3,1,6,4,7,5},
      {0,2,1,3,4,6,5,7},
      {1,0,5,4,3,2,7,6},
      {0,1,4,5,2,3,6,7},
      {7,3,5,1,6,2,4,0},
      {3,7,1,5,2,6,0,4},
      {7,6,3,2,5,4,1,0},
      {6,7,2,3,4,5,0,1},
      {2,6,0,4,3,7,1,5},
      {6,2,4,0,7,3,5,1},
      {4,6,5,7,0,2,1,3},
      {6,4,7,5,2,0,3,1},
      {7,5,6,4,3,1,2,0},
      {5,7,4,6,1,3,0,2},
      {0,4,1,5,2,6,3,7},
      {4,0,5,1,6,2,7,3},
      {3,1,7,5,2,0,6,4},
      {1,3,5,7,0,2,4,6},
      {5,7,1,3,4,6,0,2},
      {7,5,3,1,6,4,2,0},
      {3,7,2,6,1,5,0,4},
      {7,3,6,2,5,1,4,0},
      {0,2,4,6,1,3,5,7},
      {2,0,6,4,3,1,7,5},
      {5,1,4,0,7,3,6,2},
      {1,5,0,4,3,7,2,6},
      {6,2,7,3,4,0,5,1},
      {2,6,3,7,0,4,1,5},
      {6,4,2,0,7,5,3,1},
      {4,6,0,2,5,7,1,3}
      };

  constexpr int num_isometries = 48;

}//end namespace internal
}//end namespace Mesh_3
}//end namespace CGAL

#endif // CGAL_MESH_3_FEATURES_DETECTION_CUBE_ISOMETRIES_H
