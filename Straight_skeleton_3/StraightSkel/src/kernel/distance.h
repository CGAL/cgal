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
 * @file   kernel/distance.h
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#ifndef DISTANCE_H
#define DISTANCE_H

#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Plane3.h"
#include "kernel/Line3.h"

namespace kernel {

double distance(const Point2* p, const Point2* q);
double distance(const Line2* line, const Point2* point);

double distance(const Point3* p, const Point3* q);
double distance(const Plane3* plane, const Point3* point);
double distance(const Line3* line, const Point3* point);

}

#endif /* DISTANCE_H */

