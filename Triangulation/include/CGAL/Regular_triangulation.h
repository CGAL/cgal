// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
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
// Author(s)    : Samuel Hornus

#ifndef CGAL_REGULAR_TRIANGULATION_H
#define CGAL_REGULAR_TRIANGULATION_H

#include <CGAL/Triangulation.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default.h>

namespace CGAL {

template< typename RTTraits, typename TDS_ = Default >
class Regular_triangulation
: public Triangulation<RTTraits,
         typename Default::Get<TDS_, Triangulation_data_structure<
                             typename Maximal_dimension<typename RTTraits::Point_d>::type,
                             Triangulation_vertex<RTTraits>,
                             Triangulation_full_cell<RTTraits> >
                    >::type >
{
    typedef typename Maximal_dimension<typename RTTraits::Point_d>::type
                                                    Maximal_dimension_;
    typedef typename Default::Get<TDS_, Triangulation_data_structure<
                         Maximal_dimension_,
                         Triangulation_vertex<RTTraits>,
                         Triangulation_full_cell<RTTraits> >
                >::type                         TDS;
    typedef Triangulation<RTTraits, TDS>        Base;
    typedef Regular_triangulation<RTTraits, TDS_>    Self;

public:
    typedef Maximal_dimension_                  Maximal_dimension;
};

} //namespace CGAL

#endif CGAL_REGULAR_TRIANGULATION_H
