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
 * @file   kernel/projection.h
 * @author Gernot Walzl
 * @date   2013-01-29
 */

#ifndef PROJECTION_H
#define PROJECTION_H

#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Line3.h"
#include "kernel/Plane3.h"

namespace kernel {

Point2* projection(const Line2* line, const Point2* point);

Point3* projection(const Plane3* plane, const Point3* point);
Point3* projection(const Line3* line, const Point3* point);

}

#endif /* PROJECTION_H */

