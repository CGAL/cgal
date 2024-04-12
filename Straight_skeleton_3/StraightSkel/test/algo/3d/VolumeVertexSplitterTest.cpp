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

#include <cmath>
#include "algo/3d/VolumeVertexSplitter.h"
#include "algo/3d/PolyhedronBuilder.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"

BOOST_AUTO_TEST_SUITE(VolumeVertexSplitterTest)

using algo::_3d::VolumeVertexSplitter;
using algo::_3d::PolyhedronBuilder;
using data::_3d::KernelFactory;
using data::_3d::Point3SPtr;
using data::_3d::PolyhedronSPtr;

BOOST_AUTO_TEST_CASE(testCalcSurfaceArea) {
    Point3SPtr p1 = KernelFactory::createPoint3(-10.0, -10.0, -10.0);
    Point3SPtr p2 = KernelFactory::createPoint3(10.0, 10.0, -10.0);
    Point3SPtr p3 = KernelFactory::createPoint3(10.0, -10.0, 10.0);
    Point3SPtr p4 = KernelFactory::createPoint3(-10.0, 10.0, 10.0);
    PolyhedronSPtr polyhedron =
            PolyhedronBuilder::makeTetrahedron(p1, p2, p3, p4);
    double a = sqrt(20*20 + 20*20);
    double expected = sqrt(3) * a * a;
    CGAL::FT result = VolumeVertexSplitter::calcSurfaceArea(polyhedron);
    const double e = 0.001;
    BOOST_CHECK_CLOSE(expected, result, e);
}

BOOST_AUTO_TEST_SUITE_END()
