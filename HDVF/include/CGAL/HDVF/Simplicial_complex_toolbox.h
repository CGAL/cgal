// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_FILTRATION_CORE_H
#define CGAL_FILTRATION_CORE_H

#include <iostream>
#include "tools_io.hpp"
#include "Simplex.hpp"
#include "Abstract_simplicial_chain_complex.h"

namespace CGAL {
namespace HDVF {

typedef int CoefficientType;
typedef Abstract_simplicial_chain_complex<CoefficientType> ComplexType;

template <typename CoefficientType>
ComplexType& generate_Q_acyclic_complex(int n, std::set<int> subset_A)
{
    int q = subset_A.size();
}

} // end namespace HDVF
} // end namespace CGAL

#endif //CGAL_FILTRATION_CORE_H
