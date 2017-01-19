// Copyright (c) 2015 Utrecht University (The Netherlands).
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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_GET_DIMENSION_TAG_H
#define CGAL_INTERNAL_GET_DIMENSION_TAG_H

#include <CGAL/license/Spatial_searching.h>


#include <CGAL/Dimension.h>
#include <boost/mpl/has_xxx.hpp>

namespace CGAL{

namespace internal{

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_dimension_tag,Dimension,false)

  template <class T, bool has_dim = Has_dimension_tag<T>::value>
  struct Get_dimension_tag
  {
    typedef typename T::Dimension Dimension;
  };

  template <class T>
  struct Get_dimension_tag<T, false>{
    typedef Dynamic_dimension_tag Dimension;
  };

} } // end of namespace internal::CGAL

#endif //CGAL_INTERNAL_GET_DIMENSION_TAG_H
