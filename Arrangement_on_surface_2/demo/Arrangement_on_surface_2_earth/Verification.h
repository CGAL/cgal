// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef VERIFICATION_H
#define VERIFICATION_H

#include <qvector3d.h>
#include <qmatrix4x4.h>

#include "Camera.h"

// This class contains code to verify or compute certain hypotheses
//
class Verification {
public:
  // Use this to find the approximate of the true minimum projected error.
  // we are not using this complicated method, but provide it for completeness.
  static void find_minimum_projected_error_on_sphere(float we, Camera& cam,
                                                     int vp_width,
                                                     int vp_height);

  // verify that the node (180.0, -84.71338) in Antarctica is redundant
  static void verify_antarctica_node_is_redundant();
};

#endif
