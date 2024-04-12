// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   algo/2d/ptrs.h
 * @author Gernot Walzl
 * @date   2012-02-06
 */

#ifndef ALGO_2D_PTRS_H
#define ALGO_2D_PTRS_H

#include "smarter_ptr.h"

namespace algo { namespace _2d {

class SimpleStraightSkel;
class FastStraightSkel;
class SkelMeshGenerator;

typedef SHARED_PTR<SimpleStraightSkel> SimpleStraightSkelSPtr;
typedef WEAK_PTR<SimpleStraightSkel> SimpleStraightSkelWPtr;
typedef SHARED_PTR<FastStraightSkel> FastStraightSkelSPtr;
typedef WEAK_PTR<FastStraightSkel> FastStraightSkelWPtr;
typedef SHARED_PTR<SkelMeshGenerator> SkelMeshGeneratorSPtr;
typedef WEAK_PTR<SkelMeshGenerator> SkelMeshGeneratorWPtr;

} }

#endif /* ALGO_2D_PTRS_H */
