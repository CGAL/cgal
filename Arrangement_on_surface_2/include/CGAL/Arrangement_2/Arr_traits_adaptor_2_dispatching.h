// Copyright (c) 2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Eric Berberich    <eric@mpi-inf.mpg.de>

#ifndef CGAL_ARR_TRAITS_ADAPTOR_2_DISPATCHING_H
#define CGAL_ARR_TRAITS_ADAPTOR_2_DISPATCHING_H

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/has_xxx.hpp>

#include <CGAL/Arr_tags.h>

/*! \file
 * Definition of the dispatching tags and classes for the traits adaptor
 */

namespace CGAL {

//! tag to specify to use a dummy implementation
struct Arr_use_dummy_tag {};

//! tag to specify to call the corresponding traits method
struct Arr_use_traits_tag {};

namespace internal {

///////////////
// combining //
///////////////

//! struct to combine results in "or"-fashion
template < class ArrSmallerImplementationTag, class ArrLargerImplementationTag >
struct Or_traits {
  
public:
  
  //! This instance's first template parameter
  typedef ArrSmallerImplementationTag   Arr_smaller_implementation_tag;
  
  //! This instance's second template parameter
  typedef ArrLargerImplementationTag  Arr_larger_implementation_tag;
  
  
private:
  
  typedef boost::mpl::bool_< true > true_;
  typedef boost::mpl::bool_< false > false_;
  
  typedef boost::mpl::if_< 
       boost::is_same< Arr_smaller_implementation_tag, Arr_use_traits_tag >,
       true_, false_ > Smaller_traits;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_larger_implementation_tag, Arr_use_traits_tag >,
       true_, false_ > Larger_traits;

public:

  //! the result type (if one side asks for traits, then ask traits!
  //! Or vice versa: If both ask for dummy, then dummy!)
  typedef typename boost::mpl::if_<
              boost::mpl::or_< Smaller_traits, Larger_traits >,
              Arr_use_traits_tag,
              Arr_use_dummy_tag >::type type;

};


////////////////
// left-right //
////////////////

namespace Parameter_space_in_x_2 {

  // Curve-end

  template < class ArrSideTag >
  struct Curve_end {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_end< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_end< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Curve

  template < class ArrSideTag >
  struct Curve {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  // Point

  template < class ArrSideTag >
  struct Point {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Point< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
} // namespace Parameter_space_in_x_2


namespace Is_on_y_identification_2 {

  // Curve

  template < class ArrSideTag >
  struct Curve {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Point

  template < class ArrSideTag >
  struct Point {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Point< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
} // namespace Is_on_y_identification_2


namespace Compare_y_on_boundary_2 {

 // Poitns
 template < class ArrSideTag >
  struct Points {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Points< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Points< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

} // namespace Compare_y_on_boundary_2


namespace Compare_y_near_boundary_2 {

// Curve-ends

  template < class ArrSideTag >
  struct Curve_ends {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

} // Compare_y_near_boundary_2



template < class ArrLeftSideTag, class ArrRightSideTag >
struct Arr_left_right_implementation_dispatch {

public:
  
  //! This instance's first template parameter
  typedef ArrLeftSideTag   Left_side_category;
  
  //! This instance's second template parameter
  typedef ArrRightSideTag  Right_side_category;
  
public:
  
  //! tag type for Parameter_space_in_x_2 (curve-end signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Parameter_space_in_x_2::Curve_end< Left_side_category >::type,
    typename 
    internal::Parameter_space_in_x_2::Curve_end< Right_side_category >::type 
  >::type
  Parameter_space_in_x_2_curve_end_tag;

  //! tag type for Parameter_space_in_x_2 (curve signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Parameter_space_in_x_2::Curve< Left_side_category >::type,
    typename 
    internal::Parameter_space_in_x_2::Curve< Right_side_category >::type 
  >::type
  Parameter_space_in_x_2_curve_tag;
  
  //! tag type for Parameter_space_in_x_2 (point signature)
  typedef typename internal::Or_traits<
    typename
    internal::Parameter_space_in_x_2::Point< Left_side_category >::type,
    typename
    internal::Parameter_space_in_x_2::Point< Right_side_category >::type >::type
  Parameter_space_in_x_2_point_tag;

  //! tag type for Is_on_y_identification_2 (curve signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Is_on_y_identification_2::Curve< Left_side_category >::type,
    typename 
    internal::Is_on_y_identification_2::Curve< Right_side_category >::type 
  >::type
  Is_on_y_identification_2_curve_tag;

  //! tag type for Is_on_y_identification_2 (point signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Is_on_y_identification_2::Point< Left_side_category >::type,
    typename 
    internal::Is_on_y_identification_2::Point< Right_side_category >::type 
  >::type
  Is_on_y_identification_2_point_tag;

  //! tag type for Compare_y_on_boundary_2 (points signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_y_on_boundary_2::Points< Left_side_category >::type,
    typename 
    internal::Compare_y_on_boundary_2::Points< Right_side_category >::type 
  >::type
  Compare_y_on_boundary_2_points_tag;

  //! tag type for Compare_y_near_boundary_2 (curve-ends signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_y_near_boundary_2::Curve_ends< Left_side_category >::type,
    typename 
    internal::Compare_y_near_boundary_2::Curve_ends< Right_side_category >::type 
  >::type
  Compare_y_near_boundary_2_curve_ends_tag;

}; // left-right-dispatch


////////////////
// bottom-top //
////////////////

namespace Parameter_space_in_y_2 {

  // Curve-end

  template < class ArrSideTag >
  struct Curve_end {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_end< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_end< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Curve

  template < class ArrSideTag >
  struct Curve {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  // Point

  template < class ArrSideTag >
  struct Point {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Point< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
} // namespace Parameter_space_in_y_2


namespace Is_on_x_identification_2 {

  // Curve

  template < class ArrSideTag >
  struct Curve {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Point

  template < class ArrSideTag >
  struct Point {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Point< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
} // namespace Is_on_x_identification_2


namespace Compare_x_at_limit_2 {

  // Point_curve-end

  template < class ArrSideTag >
  struct Point_curve_end {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  // Curve-ends

  template < class ArrSideTag >
  struct Curve_ends {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Curve_ends< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

} // Compare_x_at_limit_2

namespace Compare_x_near_limit_2 {

  // Curve-ends

  template < class ArrSideTag >
  struct Curve_ends {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

} // Compare_x_near_limit_2


namespace Compare_x_on_boundary_2 {

 // Points
 template < class ArrSideTag >
  struct Points {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Points< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Points< Arr_contracted_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Points< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Point_curve-end

  template < class ArrSideTag >
  struct Point_curve_end {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point_curve_end< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  // Curve-ends

  template < class ArrSideTag >
  struct Curve_ends {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };

  template <>
  struct Curve_ends< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_closed_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_identified_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
} // namespace Compare_x_on_boundary_2


namespace Compare_x_near_boundary_2 {

  // Curve-ends

  template < class ArrSideTag >
  struct Curve_ends {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_oblivious_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_open_side_tag > {
    typedef Arr_use_dummy_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_ends< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

} // Compare_y_near_boundary_2



template < class ArrBottomSideTag, class ArrTopSideTag >
struct Arr_bottom_top_implementation_dispatch {

public:
  
  //! This instance's first template parameter
  typedef ArrBottomSideTag   Bottom_side_category;
  
  //! This instance's second template parameter
  typedef ArrTopSideTag  Top_side_category;
  
public:

  //! tag type for Parameter_space_in_y_2 (curve-end signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Parameter_space_in_y_2::Curve_end< Bottom_side_category >::type,
    typename 
    internal::Parameter_space_in_y_2::Curve_end< Top_side_category >::type 
  >::type
  Parameter_space_in_y_2_curve_end_tag;

  //! tag type for Parameter_space_in_y_2 (curve signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Parameter_space_in_y_2::Curve< Bottom_side_category >::type,
    typename 
    internal::Parameter_space_in_y_2::Curve< Top_side_category >::type 
  >::type
  Parameter_space_in_y_2_curve_tag;
  
  //! tag type for Parameter_space_in_y_2 (point signature)
  typedef typename internal::Or_traits<
    typename
    internal::Parameter_space_in_y_2::Point< Bottom_side_category >::type,
    typename
    internal::Parameter_space_in_y_2::Point< Top_side_category >::type >::type
  Parameter_space_in_y_2_point_tag;
  
  //! tag type for Is_on_x_identification_2 (curve signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Is_on_x_identification_2::Curve< Bottom_side_category >::type,
    typename 
    internal::Is_on_x_identification_2::Curve< Top_side_category >::type 
  >::type
  Is_on_x_identification_2_curve_tag;

  //! tag type for Is_on_x_identification_2 (point signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Is_on_x_identification_2::Point< Bottom_side_category >::type,
    typename 
    internal::Is_on_x_identification_2::Point< Top_side_category >::type 
  >::type
  Is_on_x_identification_2_point_tag;

  //! tag type for Compare_x_at_limit_2 (point-curve-end signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_at_limit_2::Point_curve_end< Bottom_side_category >
    ::type,
    typename 
    internal::Compare_x_at_limit_2::Point_curve_end< Top_side_category >::type 
  >::type
  Compare_x_at_limit_2_point_curve_end_tag;

  //! tag type for Compare_x_at_limit_2 (curve-ends signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_at_limit_2::Curve_ends< Bottom_side_category >::type,
    typename 
    internal::Compare_x_at_limit_2::Curve_ends< Top_side_category >::type 
  >::type
  Compare_x_at_limit_2_curve_ends_tag;

  //! tag type for Compare_x_near_limit_boundary_2 (curve-ends signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_near_limit_2::Curve_ends< Bottom_side_category >
    ::type,
    typename 
    internal::Compare_x_near_limit_2::Curve_ends< Top_side_category >::type 
  >::type
  Compare_x_near_limit_2_curve_ends_tag;

  //! tag type for Compare_x_on_boundary_2 (points signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_on_boundary_2::Points< Bottom_side_category >::type,
    typename 
    internal::Compare_x_on_boundary_2::Points< Top_side_category >::type 
  >::type
  Compare_x_on_boundary_2_points_tag;

  //! tag type for Compare_x_on_boundary_2 (point_curve-end signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_on_boundary_2::Point_curve_end< Bottom_side_category >::type,
    typename 
    internal::Compare_x_on_boundary_2::Point_curve_end< Top_side_category >::type 
  >::type
  Compare_x_on_boundary_2_point_curve_end_tag;

  //! tag type for Compare_x_on_boundary_2 (curve-ends signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_on_boundary_2::Curve_ends< Bottom_side_category >::type,
    typename 
    internal::Compare_x_on_boundary_2::Curve_ends< Top_side_category >::type 
  >::type
  Compare_x_on_boundary_2_curve_ends_tag;
  
  //! tag type for Compare_x_near_boundary_2 (curve-ends signature)
  typedef typename internal::Or_traits<
    typename 
    internal::Compare_x_near_boundary_2::Curve_ends< Bottom_side_category >::type,
    typename 
    internal::Compare_x_near_boundary_2::Curve_ends< Top_side_category >::type 
  >::type
  Compare_x_near_boundary_2_curve_ends_tag;

  }; // bottom-top-dispatch

} // namespace internal

} //namespace CGAL

#endif // CGAL_ARR_TRAITS_ADAPTOR_2_DISPATCHING_H
