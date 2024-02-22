// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#pragma once

#include "geometry_basics.h"
#include "certificate.h"
#include "curve.h"

#include <array>

class FrechetAbstract
{
public:
	virtual ~FrechetAbstract() {}
	virtual bool lessThan(distance_t distance, Curve const& curve1, Curve const& curve2) = 0;
	virtual Certificate&  computeCertificate() = 0;

	// yes, this is ugly...
	virtual void setRules(std::array<bool,5> const& enable) {}
	virtual void setPruningLevel(int pruning_level) {};
};
