// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H
#define CGAL_ARR_CONSOLIDATED_CURVE_DATA_TRAITS_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * Definition of the Arr_consolidated_curve_data_traits_2<Traits,Data> class.
 */

#include<CGAL/Arr_curve_data_traits_2.h>
#include<CGAL/Arr_geometry_traits/Consolidated_curve_data_aux.h>

namespace CGAL {

/*! \class
 * A generic traits class for maintaining an arrangement of curves that have
 * an extra data field. This traits class is templated with a Data class an
 * an ordinary traits class which is also used as a based traits class to
 * inherit from. It extracts the original Curve_2 and X_monotone_curve_2 types
 * from the ordinary traits class, and redefines them to have Data as an extra
 * field in the Curve_2 type, and a container of Data objects for the extended
 * X_monotone_curve_2 type.
 * The Data field is updated when the curves are converted from Curve_2 to
 * X_monotone_curve_2, and when the X_monotone_curve_2 curves are split.
 * When two x-monotone curves overlap, their data containers are consolidated
 * and attached to the resulting subcurve.
 * All other functors are inherited from the base ordinary traits class.
 */
template <class Traits_, class Data_>
class Arr_consolidated_curve_data_traits_2 :
  public Arr_curve_data_traits_2<Traits_,
                                 _Unique_list<Data_>, 
                                 _Consolidate_unique_lists<Data_>,
                                 Data_>
{
private:

  typedef Arr_curve_data_traits_2<Traits_,
                                  _Unique_list<Data_>, 
                                  _Consolidate_unique_lists<Data_>,
                                  Data_>            Base;

public:

  typedef Traits_                                     Base_traits_2;
  typedef Data_                                       Data;
  typedef _Unique_list<Data_>                         Data_container;
  typedef typename Data_container::const_iterator     Data_iterator;
  typedef typename Data_container::const_iterator     Data_const_iterator;

  typedef typename Base::Curve_2                      Curve_2;
  typedef typename Base_traits_2::Curve_2             Base_curve_2;
  typedef typename Base::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Base_traits_2::X_monotone_curve_2  Base_x_monotone_curve_2;
  typedef typename Base_traits_2::Point_2             Point_2;
  typedef typename Base_traits_2::Multiplicity        Multiplicity;

  typedef typename Base_traits_2::Has_left_category   Has_left_category;
  
  typedef typename Base_traits_2::Has_merge_category  Base_has_merge_category;
  typedef Tag_true                                    Has_merge_category;
  typedef typename Base_traits_2::Has_do_intersect_category
                                                      Has_do_intersect_category;

  // Base_traits_2 is Arr_curve_data_traits that already completes
  // incomplete tags
  typedef typename Base_traits_2::Left_side_category
                                                      Left_side_category;
  typedef typename Base_traits_2::Bottom_side_category
                                                      Bottom_side_category;
  typedef typename Base_traits_2::Top_side_category
                                                      Top_side_category;
  typedef typename Base_traits_2::Right_side_category
                                                      Right_side_category;
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
