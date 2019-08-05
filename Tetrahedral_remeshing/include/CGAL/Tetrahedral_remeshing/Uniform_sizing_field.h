// Copyright (c) 2019 GeometryFactory (France).
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
//
// Author(s)     : Jane Tournois

#ifndef CGAL_UNIFORM_SIZING_FIELD_H
#define CGAL_UNIFORM_SIZING_FIELD_H

#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>

namespace CGAL
{
  template <class Kernel>
  class Uniform_sizing_field : Sizing_field<Kernel>
  {
  private:
    typedef Sizing_field<Kernel>        Base;
  public:
    typedef typename Base::FT         FT;
    typedef typename Base::Point_3    Point_3;

    Uniform_sizing_field(const FT& size)
      : m_size(size)
    {}

    FT operator()(const Point_3&) const
    {
      return m_size;
    }
  private:
    FT m_size;
  };

}//end namespace CGAL

#endif //CGAL_UNIFORM_SIZING_FIELD_H
