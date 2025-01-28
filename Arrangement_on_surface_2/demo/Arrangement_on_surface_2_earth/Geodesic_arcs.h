// Copyright (c) 2023, 2024 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Engin Deniz Diktas <denizdiktas@gmail.com>

#ifndef GEODESIC_ARCS_H
#define GEODESIC_ARCS_H

#include <vector>
#include <qvector3d.h>

#include "Kml_reader.h"

class Geodesic_arcs {
public:
  using Approx_arcs = std::vector<std::vector<QVector3D>>;

  Approx_arcs get_approx_arcs(double error);

  // generate approximate arcs from KML data
  Approx_arcs get_approx_arcs(const Kml::Placemark& placemark, double error);
  Approx_arcs get_approx_arcs(const Kml::Placemarks& placemarks, double error);
};


#endif
