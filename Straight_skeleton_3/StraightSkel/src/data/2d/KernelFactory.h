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
 * @file   data/2d/KernelFactory.h
 * @author Gernot Walzl
 * @date   2011-11-25
 */

#ifndef DATA_2D_KERNELFACTORY_H
#define DATA_2D_KERNELFACTORY_H

#include "config.h"
#include "cgal_kernel.h"

#ifndef USE_CGAL
    #include "kernel/Point2.h"
    #include "kernel/Line2.h"
    #include "kernel/Segment2.h"
    #include "kernel/Vector2.h"
#endif

#include "debug.h"
#include "data/2d/ptrs.h"

namespace data { namespace _2d {

class KernelFactory {
public:
    virtual ~KernelFactory();

    static Point2SPtr createPoint2(CGAL::FT x, CGAL::FT y);
    static Point2SPtr createPoint2(const Point2& point);

    static Segment2SPtr createSegment2(Point2SPtr src, Point2SPtr dst);
    static Segment2SPtr createSegment2(const Segment2& seg);

    static Line2SPtr createLine2(Point2SPtr p, Point2SPtr q);
    static Line2SPtr createLine2(const Line2& line);
    static Line2SPtr createLine2(Segment2SPtr seg);
    static Line2SPtr createLine2(Point2SPtr p, Vector2SPtr direction);

    static Vector2SPtr createVector2(CGAL::FT x, CGAL::FT y);
    static Vector2SPtr createVector2(const Vector2& vector);
    static Vector2SPtr createVector2(Point2SPtr point);
    static Vector2SPtr createVector2(Line2SPtr line);

protected:
    KernelFactory();
};

} }

#endif /* DATA_2D_KERNELFACTORY_H */

