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
 * @file   kernel/Ray2.h
 * @author Gernot Walzl
 * @date   2012-02-13
 */

#ifndef RAY2_H
#define RAY2_H

#include "kernel/Point2.h"
#include "kernel/Vector2.h"

namespace kernel {

class Ray2 {
public:
    Ray2();
    Ray2(const Point2& point, const Vector2& direction);
    Ray2(const Ray2& orig);
    virtual ~Ray2();
protected:
    Point2* point_;
    Vector2* direction_;
};

}

#endif /* RAY2_H */

