// Copyright (c) 2007, 2009 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Eric Berberich    <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_TAGS_H
#define CGAL_ARR_TAGS_H

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>

/*! \file
 * Definition of the tags for the arrangement package.
 */

CGAL_BEGIN_NAMESPACE

struct Arr_no_boundary_tag {};
struct Arr_has_boundary_tag {};
struct Arr_bounded_boundary_tag : public virtual Arr_has_boundary_tag {};
struct Arr_unbounded_boundary_tag : public virtual Arr_has_boundary_tag {};
struct Arr_all_boundary_tag : public virtual Arr_bounded_boundary_tag,
                              public virtual Arr_unbounded_boundary_tag {};

struct Arr_boundary_side_tag {};
struct Arr_oblivious_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_open_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_closed_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_contracted_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_identified_side_tag : public virtual Arr_boundary_side_tag {};

struct Arr_all_sides_oblivious_tag :
  public boost::mpl::bool_< true > {};
struct Arr_not_all_sides_oblivious_tag :
  public boost::mpl::bool_< false > {};

struct Arr_all_sides_non_open_tag : 
  public virtual boost::mpl::bool_< true > {};
struct Arr_not_all_sides_non_open_tag : 
  public virtual boost::mpl::bool_< false > {};


/*!\brief Struct to determine whether all side tags are "oblivious"
 */
template < class ArrLeftSideTag, class ArrBottomSideTag, 
           class ArrTopSideTag, class ArrRightSideTag >
struct Arr_are_all_sides_oblivious_tag {

public:

  //! This instance's first template parameter
  typedef ArrLeftSideTag   Arr_left_side_tag;
  
  //! This instance's second template parameter
  typedef ArrBottomSideTag Arr_bottom_side_tag;
  
  //! This instance's third template parameter
  typedef ArrTopSideTag    Arr_top_side_tag;
  
  //! This instance's fourth template parameter
  typedef ArrRightSideTag  Arr_right_side_tag;
  
private:
  
  typedef boost::mpl::bool_< true > true_;
  typedef boost::mpl::bool_< false > false_;
  
  typedef boost::mpl::if_< 
       boost::is_same< Arr_left_side_tag, Arr_oblivious_side_tag >,
       true_, false_ > 
  Left_oblivious;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_bottom_side_tag, Arr_oblivious_side_tag >,
       true_, false_ > 
  Bottom_oblivious;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_top_side_tag, Arr_oblivious_side_tag >,
       true_, false_ > 
  Top_oblivious;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_right_side_tag, Arr_oblivious_side_tag >,
       true_, false_ > 
  Right_oblivious;
  
public:
  
  /*!\brief
   * boolean tag that is Arr_all_sides_oblivious_tag if all sides are 
   * oblivious, otherwise Arr_not_all_sides_oblivious_tag
   */
  typedef typename boost::mpl::if_< 
                           boost::mpl::and_< Left_oblivious, Bottom_oblivious, 
                                             Top_oblivious, Right_oblivious >,
                           Arr_all_sides_oblivious_tag,
                           Arr_not_all_sides_oblivious_tag >::type result;

};

/*!\brief Struct to determine whether all side tags are "non-open"
 */
template < class ArrLeftSideTag, class ArrBottomSideTag, 
           class ArrTopSideTag, class ArrRightSideTag >
struct Arr_are_all_sides_non_open_tag {

public:

  //! This instance's first template parameter
  typedef ArrLeftSideTag   Arr_left_side_tag;
  
  //! This instance's second template parameter
  typedef ArrBottomSideTag Arr_bottom_side_tag;
  
  //! This instance's third template parameter
  typedef ArrTopSideTag    Arr_top_side_tag;
  
  //! This instance's fourth template parameter
  typedef ArrRightSideTag  Arr_right_side_tag;
  
private:
  
  typedef boost::mpl::bool_< true > true_;
  typedef boost::mpl::bool_< false > false_;
  
  typedef boost::mpl::if_< 
       boost::is_same< Arr_left_side_tag, Arr_open_side_tag >,
       true_, false_ > 
  Left_open;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_bottom_side_tag, Arr_open_side_tag >,
       true_, false_ > 
  Bottom_open;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_top_side_tag, Arr_open_side_tag >,
       true_, false_ > 
  Top_open;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_right_side_tag, Arr_open_side_tag >,
       true_, false_ > 
  Right_open;
  
public:
  
  /*!\brief
   * boolean tag that is Arr_all_sides_non_open_tag if all sides are non-open, 
   * otherwise Arr_not_all_sides_non_open_tag
   */
  typedef typename boost::mpl::if_<
      boost::mpl::and_< boost::mpl::not_< Left_open >, 
                        boost::mpl::not_< Bottom_open >, 
                        boost::mpl::not_< Top_open >, 
                        boost::mpl::not_< Right_open > >,
      Arr_all_sides_non_open_tag,
      Arr_not_all_sides_non_open_tag >::type result;
};


/*!\brief Struct to check consistent tagging of identifications
 */
template < class ArrLeftSideTag, class ArrBottomSideTag, 
           class ArrTopSideTag, class ArrRightSideTag >
struct Arr_sane_identified_tagging {

public:
  
  //! This instance's first template parameter
  typedef ArrLeftSideTag   Arr_left_side_tag;
  
  //! This instance's second template parameter
  typedef ArrBottomSideTag Arr_bottom_side_tag;
  
  //! This instance's third template parameter
  typedef ArrTopSideTag    Arr_top_side_tag;
  
  //! This instance's fourth template parameter
  typedef ArrRightSideTag  Arr_right_side_tag;
  
private:
  
  typedef boost::mpl::bool_< true > true_;
  typedef boost::mpl::bool_< false > false_;
  
  typedef boost::mpl::if_< 
       boost::is_same< Arr_left_side_tag, Arr_identified_side_tag >,
       true_, false_ > 
  Left_identified;

  typedef boost::mpl::if_<
       boost::is_same< Arr_bottom_side_tag, Arr_identified_side_tag >,
       true_, false_ > 
  Bottom_identified;

  typedef boost::mpl::if_< 
       boost::is_same< Arr_top_side_tag, Arr_identified_side_tag >,
       true_, false_ > 
  Top_identified;
  
  typedef boost::mpl::if_< 
       boost::is_same< Arr_right_side_tag, Arr_identified_side_tag >,
       true_, false_ > 
  Right_identified;

  typedef boost::mpl::and_< Left_identified, Right_identified > LR_identified;

  typedef boost::mpl::and_< Bottom_identified, Top_identified > BT_identified;

  typedef boost::mpl::and_< boost::mpl::not_< Left_identified>,
                            boost::mpl::not_< Right_identified > >
  LR_non_identified;

  typedef boost::mpl::and_< boost::mpl::not_< Bottom_identified >, 
                            boost::mpl::not_< Top_identified > > 
  BT_non_identified;
  
  typedef boost::mpl::or_< LR_identified, LR_non_identified > LR_ok;
  typedef boost::mpl::or_< BT_identified, BT_non_identified > BT_ok;
  
public:
  
  /*!\brief
   * boolean tag that is bool_<true> if opposite sides are either 
   * both identified or both non-identified,
   * otherwise bool_<false>
   */
  typedef boost::mpl::and_< LR_ok, BT_ok > result;

};


/*!\brief Struct to determine whether the tags in a geometry traits
 * and a topology traits match.
 */
template < class GeometryTraits_2, class TopologyTraits_2 >
struct Arr_same_tags_in_traits_classes {

public:

  //! this instance's first template parameter
  typedef GeometryTraits_2 Geometry_traits_2;

  //! this instance's second template parameter
  typedef TopologyTraits_2 Topology_traits_2;

private:

  typedef typename Geometry_traits_2::Arr_left_side_tag   Geo_left_side_tag;
  typedef typename Geometry_traits_2::Arr_bottom_side_tag Geo_bottom_side_tag;
  typedef typename Geometry_traits_2::Arr_top_side_tag    Geo_top_side_tag;
  typedef typename Geometry_traits_2::Arr_right_side_tag  Geo_right_side_tag;

  typedef typename Topology_traits_2::Arr_left_side_tag   Top_left_side_tag;
  typedef typename Topology_traits_2::Arr_bottom_side_tag Top_bottom_side_tag;
  typedef typename Topology_traits_2::Arr_top_side_tag    Top_top_side_tag;
  typedef typename Topology_traits_2::Arr_right_side_tag  Top_right_side_tag;
  
  typedef boost::mpl::bool_< true > true_;
  typedef boost::mpl::bool_< false > false_;
  
  typedef boost::mpl::if_< 
       boost::is_same< Geo_left_side_tag, Top_left_side_tag >,
       true_, false_ > 
  Left_equal;

 typedef boost::mpl::if_< 
       boost::is_same< Geo_bottom_side_tag, Top_bottom_side_tag >,
       true_, false_ > 
  Bottom_equal;

 typedef boost::mpl::if_< 
       boost::is_same< Geo_top_side_tag, Top_top_side_tag >,
       true_, false_ > 
  Top_equal;

 typedef boost::mpl::if_< 
       boost::is_same< Geo_right_side_tag, Top_right_side_tag >,
       true_, false_ > 
  Right_equal;

public:
  
  /*!\brief
   * boolean tag that is bool_<true> if side tags in given traits classes 
   * match, otherwise bool_<false>
   */
  typedef boost::mpl::and_< Left_equal, Bottom_equal, 
                            Top_equal, Right_equal > result;
};


struct Arr_use_implementation_tag {};
struct Arr_use_default_tag : public virtual Arr_use_implementation_tag {};
struct Arr_use_traits_tag  : public virtual Arr_use_implementation_tag {};


namespace CGALi {

namespace Parameter_space_in_x_2 {

  // Curve-end

  template < class ArrSideTag >
  struct Curve_end {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve_end< Arr_oblivious_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve_end< Arr_open_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve_end< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Curve

  template < class ArrSideTag >
  struct Curve {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve< Arr_oblivious_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve< Arr_open_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve< Arr_contracted_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Curve< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Curve< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };

  // Point

  template < class ArrSideTag >
  struct Point {
    typedef Arr_use_default_tag type;
  };

  template <>
  struct Point< Arr_oblivious_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Point< Arr_open_side_tag > {
    typedef Arr_use_default_tag type;
  };
  
  template <>
  struct Point< Arr_contracted_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_closed_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  
  template <>
  struct Point< Arr_identified_side_tag > {
    typedef Arr_use_traits_tag type;
  };
  

} // namespace Parameter_space_in_x_2

// TODO missing functors + signature 
// TODO move to another file

} // namespace CGALi

template < class ArrLeftSideTag, class ArrRightSideTag >
struct Arr_left_right_implementation_dispatch {

public:
  
  //! This instance's first template parameter
  typedef ArrLeftSideTag   Arr_left_side_tag;
  
  //! This instance's second template parameter
  typedef ArrRightSideTag  Arr_right_side_tag;
  
private:

  //! struct to combine results in "or"-fashion
  template < class ArrLeftImplementationTag, class ArrRightImplementationTag >
  struct Or_traits {
    
  public:
    
    //! This instance's first template parameter
    typedef ArrLeftImplementationTag   Arr_left_implementation_tag;
    
    //! This instance's second template parameter
    typedef ArrRightImplementationTag  Arr_right_implementation_tag;
    
  public:
    
    //! the result type (if one side asks for traits, then ask traits!
    //! Or vice versa: If both ask for default, then default!)
    typedef boost::mpl::if_<
//              boost::mpl::or_< 
//                boost::is_same< Arr_left_implementation_tag, 
//                                Arr_use_traits_tag >,
//                boost::is_same< Arr_right_implementation_tag, 
//                                Arr_use_traits_tag >
//              >,
              boost::mpl::bool_<true>,
              Arr_use_traits_tag,
              Arr_use_default_tag >::type type;
  };
  
public:
  
  //! tag type for Parameter_space_in_x_2 (curve-end signature)
  typedef typename Or_traits<
    typename 
    CGALi::Parameter_space_in_x_2::Curve_end< Arr_left_side_tag >::type,
    typename 
    CGALi::Parameter_space_in_x_2::Curve_end< Arr_right_side_tag >::type 
  >::type
  Parameter_space_in_x_2_curve_end_tag;

  //! tag type for Parameter_space_in_x_2 (curve signature)
  typedef typename Or_traits<
    typename 
    CGALi::Parameter_space_in_x_2::Curve< Arr_left_side_tag >::type,
    typename 
    CGALi::Parameter_space_in_x_2::Curve< Arr_right_side_tag >::type 
  >::type
  Parameter_space_in_x_2_curve_tag;
  
  //! tag type for Parameter_space_in_x_2 (point signature)
  typedef typename Or_traits<
    typename
    CGALi::Parameter_space_in_x_2::Point< Arr_left_side_tag >::type,
    typename
    CGALi::Parameter_space_in_x_2::Point< Arr_right_side_tag >::type >::type
  Parameter_space_in_x_2_point_tag;
  
  // TODO missing functors + signatures

};

CGAL_END_NAMESPACE

#endif

