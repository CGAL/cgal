//
// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s): Saurabh Singh <ssingh@cs.iitr.ac.in>

#ifndef ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
#define ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H

#include <vector>
#include <string>
#include <complex>
#include <boost/optional.hpp>
#include <CGAL/Polynomial_traits_d.h>

template <typename Polynomial_2>
class AlgebraicCurveParser {
	using Traits = CGAL::Polynomial_traits_d<Polynomial_2>;
	using Coefficient = typename Traits::Innermost_coefficient_type;

public:
    explicit AlgebraicCurveParser(std::string expression);
    bool validateExpression();
	boost::optional<Polynomial_2> parse();

private:
    std::string expression;
};

#endif //ARRANGEMENT_ON_SURFACE_2_DEMO_ALGEBRAICCURVEPARSERNEW_H
