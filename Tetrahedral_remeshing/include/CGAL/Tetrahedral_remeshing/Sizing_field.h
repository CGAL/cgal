// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
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

#ifndef CGAL_SIZING_FIELD_H
#define CGAL_SIZING_FIELD_H

#include <CGAL/license/Tetrahedral_remeshing.h>

namespace CGAL
{
/*!
* Sizing field virtual class
*/
template <class Kernel>
class Sizing_field
{
public:
  typedef Kernel                      K;
  typedef typename Kernel::FT         FT;
  typedef typename Kernel::Point_3    Point_3;

public:
  virtual FT operator()(const Point_3& p) const = 0;
};

}//end namespace CGAL

#endif //CGAL_SIZING_FIELD_H
