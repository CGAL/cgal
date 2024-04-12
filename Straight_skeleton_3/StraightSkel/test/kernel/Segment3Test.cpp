// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#include <boost/test/unit_test.hpp>

#include "kernel/Segment3.h"

using kernel::Segment3;
using kernel::Point3;
using kernel::Line3;

BOOST_AUTO_TEST_SUITE(Segment3Test)

BOOST_AUTO_TEST_CASE(testConstructor) {
    Point3 p(1.0, 2.0, 3.0);
    Point3 q(3.0, 3.0, 3.0);

    Segment3 s(p,q);

    BOOST_CHECK_EQUAL(&p, &s.getP());
    BOOST_CHECK_EQUAL(&q, &s.getQ());
}

BOOST_AUTO_TEST_SUITE_END()

