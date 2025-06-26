// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2016 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#ifndef CGAL_HAS_NESTED_TYPE_POINT_H
#define CGAL_HAS_NESTED_TYPE_POINT_H

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/identity.hpp>

namespace CGAL {

  namespace internal {

    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Point, Point, false)

    template<typename Gt>
    struct Point_type
    {
    public:
      typedef typename Gt::Point type;
    };

    template<typename Gt, typename Other = void>
    struct Get_nested_type_Point
    {
      typedef typename boost::mpl::eval_if_c<
        Has_nested_type_Point<Gt>::value,
        Point_type<Gt>,
        boost::mpl::identity<Other>
      >::type type;
    };

  } // end namespace internal
} // end namespace CGAL

#endif // CGAL_HAS_NESTED_TYPE_POINT_H
