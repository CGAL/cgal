// Copyright (c) 2026  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : François Protais

#ifndef CGAL_MESH_SMOOTHING_3_INTERNAL_PREDICATES_PREDICATES_H
#define CGAL_MESH_SMOOTHING_3_INTERNAL_PREDICATES_PREDICATES_H

#include <CGAL/license/Mesh_smoothing_3.h>

#include <Eigen/Eigen>
#include <array>

namespace CGAL{

namespace Mesh_smoothing_3_internal {


namespace exact_predicates {

double orient3d(double pa[3], double pb[3], double pc[3], double pd[3]);

inline bool positive_tetrahedra(std::array<Eigen::Vector3d, 4> const &tetrahedra) {
    double pa[3] = {tetrahedra[0][0], tetrahedra[0][1], tetrahedra[0][2]};
    double pb[3] = {tetrahedra[1][0], tetrahedra[1][1], tetrahedra[1][2]};
    double pc[3] = {tetrahedra[2][0], tetrahedra[2][1], tetrahedra[2][2]};
    double pd[3] = {tetrahedra[3][0], tetrahedra[3][1], tetrahedra[3][2]};
    return orient3d(pa,pb,pc,pd) > 0;
}


#include "predicates_shewchuk.h"

}
} } // end of CGAL::Mesh_smoothing_3_internal namespace

#endif // CGAL_MESH_SMOOTHING_3_INTERNAL_PREDICATES_PREDICATES_H
