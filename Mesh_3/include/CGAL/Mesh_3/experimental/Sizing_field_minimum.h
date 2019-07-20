// Copyright (c) 2016 GeometryFactory (France).
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
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef CGAL_MESH_3_SIZING_FIELD_MINIMUM
#define CGAL_MESH_3_SIZING_FIELD_MINIMUM

#include <CGAL/license/Mesh_3.h>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL
{
  template <typename SizingField1, typename SizingField2>
  class Sizing_field_minimum
  {
  public:
    typedef typename SizingField1::FT       FT;
    typedef typename SizingField1::Point_3  Point_3;
    typedef typename SizingField1::Index    Index;

    BOOST_STATIC_ASSERT_MSG((
      boost::is_same<typename SizingField1::FT,
                     typename SizingField2::FT>::value),
      "FT type should be the same for both sizing fields");
    BOOST_STATIC_ASSERT_MSG((
      boost::is_same<typename SizingField1::Point_3,
                     typename SizingField2::Point_3>::value),
      "Point_3 type should be the same for both sizing fields");
    BOOST_STATIC_ASSERT_MSG((
      boost::is_same<typename SizingField1::Index,
                     typename SizingField2::Index>::value),
      "Index type should be the same for both sizing fields");

  private:
    const SizingField1* mp_size1;
    const SizingField2* mp_size2;

  public:
    Sizing_field_minimum(const SizingField1* sf1,
                         const SizingField2* sf2)
      : mp_size1(sf1)
      , mp_size2(sf2)
    {}

    FT operator()(const Point_3& p, const int dim, const Index& index) const
    {
      if(mp_size2 == 0) return (*mp_size1)(p, dim, index);
      if(mp_size1 == 0) return (*mp_size2)(p, dim, index);
      FT s1 = (*mp_size1)(p, dim, index);
      FT s2 = (*mp_size2)(p, dim, index);
      if (s1 == 0)       return s2;
      else if (s2 == 0)  return s1;
      else               return (std::min)(s1, s2);
    }

  };//Sizing_field_minimum

}//namespace CGAL

#endif // CGAL_MESH_3_SIZING_FIELD_MINIMUM
