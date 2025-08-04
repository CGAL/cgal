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
//                  Kevin Fedyna <fedyna.kevin@gmail.com>


#ifndef CGAL_OSM__BASE_H
#define CGAL_OSM__BASE_H

#include <CGAL/license/HDVF.h>

#include <ostream>

namespace CGAL {
namespace OSM {

/** \brief Chain type flag for column chain. */
const int COLUMN = 0b01;
/** \brief Chain type flag for row chain. */
const int ROW    = 0b10;

/** \brief The default type for signed integers. */
typedef int ZCoefficient;

// Class Sparse_matrix
template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
class Sparse_matrix;

// Class Sparse_chain

template <typename _CoefficientType = OSM::ZCoefficient, int _ChainTypeFlag = OSM::COLUMN>
class Sparse_chain;

// Class Bitboard

class Bitboard;

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM__BASE_H
