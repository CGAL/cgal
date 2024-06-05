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

#ifndef CGAL_UNIFORM_SIZING_FIELD_H
#define CGAL_UNIFORM_SIZING_FIELD_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing_sizing_field.h>

namespace CGAL
{
/**
 * @ingroup PkgTetrahedralRemeshingSizing
 *
 * @class Uniform_sizing_field
 * @tparam GT the geometric traits class
 *
 * The uniform (i.e. constant) sizing field for tetrahedral remeshing
 *
 * \cgalModels{RemeshingSizingField_3}
 */
template <class GT>
class Uniform_sizing_field
  : public Tetrahedral_remeshing_sizing_field<GT>
{
public:
  typedef typename GT::FT         FT;
  typedef typename GT::Point_3    Point_3;

  /** Constructor
  * @param size the target edge length for remeshing
  */
  Uniform_sizing_field(const FT& size)
    : m_size(size)
  {}

  template<typename Index>
  FT operator()(const Point_3&, const int, const Index&) const
  {
    return m_size;
  }

private:
  const FT m_size;
};

}//end namespace CGAL

#endif //CGAL_UNIFORM_SIZING_FIELD_H
