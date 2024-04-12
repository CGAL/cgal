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
 * @file   kernel/Segment2.h
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#ifndef SEGMENT2_H
#define SEGMENT2_H

#include "kernel/Point2.h"
#include "kernel/Line2.h"

namespace kernel {

class Segment2 {
public:
    Segment2(const Point2& p, const Point2& q);
    Segment2(const Segment2& orig);
    virtual ~Segment2();
    const Point2& getP() const;
    const Point2& getQ() const;
    void setP(const Point2& p);
    void setQ(const Point2& q);
    Line2 line() const;
protected:
    const Point2* p_;
    const Point2* q_;
};

}

#endif /* SEGMENT2_H */

