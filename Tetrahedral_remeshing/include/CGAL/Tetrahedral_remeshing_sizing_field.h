// Copyright (c) 2024 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_TETRAHEDRAL_REMESHING_SIZING_FIELD_H
#define CGAL_TETRAHEDRAL_REMESHING_SIZING_FIELD_H

#include <CGAL/license/Tetrahedral_remeshing.h>

namespace CGAL
{

/*!
* @ingroup PkgTetrahedralRemeshingSizing
*
* Sizing field virtual class, designed for tetrahedral remeshing
* @tparam GT the geometric traits class
*/
template <class GT>
class Tetrahedral_remeshing_sizing_field
{
public:
  typedef GT                      K;
  typedef typename GT::FT         FT;
  typedef typename GT::Point_3    Point_3;

public:
  template<typename Index>
  FT operator()(const Point_3& p, const int dim, const Index& i) const;
};

}//end namespace CGAL

#endif //CGAL_TETRAHEDRAL_REMESHING_SIZING_FIELD_H
