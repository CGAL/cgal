// Copyright(c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef AOS_TRIANGULATOR_H
#define AOS_TRIANGULATOR_H

#include <map>
#include <vector>

#include <qvector3d.h>

#include "Aos.h"

using Country_triangles_map = std::map<std::string, std::vector<QVector3D>>;

class Aos_triangulator {
public:
  static std::vector<QVector3D> get_all(Aos::Arr_handle arrh);

  static Country_triangles_map get_by_country(Aos::Arr_handle arrh,
                                              float error,
                                              std::size_t num_uniform_points);
};

#endif
