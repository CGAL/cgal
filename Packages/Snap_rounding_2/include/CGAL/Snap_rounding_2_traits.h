// ======================================================================
//
// Copyright (c) The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-24 $
// release_date  : $CGAL_Date: 2000/12/29 $
//
// file          : include/CGAL/Snap_rounding_2_traits.h
// package       : arr (1.73)
// maintainer    : Eli Packer <elip@post.tau.ac.il>
// author(s)     : Eli Packer
// coordinator   : Tel-Aviv University (Dan Halperin <danha@post.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_SR_2_TRAITS_H
#define CGAL_SR_2_TRAITS_H

#include <CGAL/basic.h>
#include <map>

#include <CGAL/leda_integer.h> 
#include "../../include/CGAL/Snap_rounding_2_utility.h"

CGAL_BEGIN_NAMESPACE

template<class base_rep>
class Snap_rounding_traits : public base_rep {

typedef typename base_rep::FT                    NT;
typedef typename base_rep::Point_2               Point_2;
typedef typename base_rep::Segment_2             Segment_2;
typedef typename base_rep::Iso_rectangle_2       Iso_rectangle_2;

public:

Snap_rounding_traits()
  {
    init_angle_appr();
  }

/*! Functor
 */
class Snap_2 {
 public:
  void operator()(Point_2 p,NT pixel_size,NT &x,NT &y)
  {
    NT x_tmp = p.x() / pixel_size;
    NT y_tmp = p.y() / pixel_size;

    x = floor(x_tmp.to_double()) * pixel_size + pixel_size / 2.0;
    y = floor(y_tmp.to_double()) * pixel_size + pixel_size / 2.0;
  }
};

Snap_2 snap_2_object() const { return Snap_2(); }

/*! Functor
 */
class Integer_grid_point_2 {
 public:
  Point_2 operator()(Point_2 p,NT pixel_size)
  {
    NT x = (p.x() - pixel_size / 2.0) / pixel_size;
    NT y = (p.y() - pixel_size / 2.0) / pixel_size;
    Point_2 out_p(x,y);

    return(out_p);
  }
};

Integer_grid_point_2 integer_grid_point_2_object() const
    { return Integer_grid_point_2(); }

/*! Functor
 */
class Segment_direction_2 {
 public:
  double operator()(Segment_2 s)
  {
    double x1 = s.source().x().to_double();
    double y1 = s.source().y().to_double();
    double x2 = s.target().x().to_double();
    double y2 = s.target().y().to_double();

    return(atan((y2 - y1)/(x2 - x1)));
  }
};

Segment_direction_2 segment_direction_2_object() const
    {return Segment_direction_2(); }

/*! Functor
 */
/*class Rotate_point_2 {
 public:
  void operator()(Point_2& p,NT angle)
  {
    int tranc_angle = int(angle.to_double() * rad_to_deg);
    NT cosine_val = angle_to_sines_appr[90 - tranc_angle],
       sine_val = angle_to_sines_appr[tranc_angle];

    Transformation_2 rotate(ROTATION, sine_val, cosine_val);

    p = rotate(p);
  }
};

Rotate_point_2 rotate_point_2_object() const {return Rotate_point_2(); }
*/
/*! Functor
 */

class Bounding_box_of_minkowski_sum_2 {
 private:
  const Snap_rounding_traits<base_rep>* _gt;
  Bounding_box_of_minkowski_sum_2(const Snap_rounding_traits<base_rep>* gt) : _gt(gt) {}

  Point_2 small_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_x_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_x_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

  Point_2 small_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p1);
    else
      return(p2);
  }

  Point_2 big_y_point(Point_2 p1,Point_2 p2)
  {
    Comparison_result c = _gt->compare_y_2_object()(p1,p2);
    if(c == SMALLER)
      return(p2);
    else
      return(p1);
  }

 public:
  Iso_rectangle_2 operator()(Segment_2 s,
                             NT unit_squere,
		             NT angle)
  {
    Point_2 ms1,ms2,ms3,ms4,ms5,ms6;// minkowski sum points

    Comparison_result cx =  _gt->compare_x_2_object()(s.source(),s.target());
    NT x1 = s.source().x(),y1 = s.source().y(),x2 =
       s.target().x(),y2 = s.target().y();

    if(cx == SMALLER) {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point_2(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point_2(x1 - 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms3 = Point_2(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms4 = Point_2(x2 + 0.6 * unit_squere,y2 - 0.6 * unit_squere);
      ms5 = Point_2(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point_2(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
    } else {
      // we use unit_squere instead of unit_squere / 2 in order to
      // find tangency points which are not supported by kd-tree
      ms1 = Point_2(x1 + 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms2 = Point_2(x1 - 0.6 * unit_squere,y1 - 0.6 * unit_squere);
      ms3 = Point_2(x1 + 0.6 * unit_squere,y1 + 0.6 * unit_squere);
      ms4 = Point_2(x2 + 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms5 = Point_2(x2 - 0.6 * unit_squere,y2 + 0.6 * unit_squere);
      ms6 = Point_2(x2 - 0.6 * unit_squere,y2 - 0.6 * unit_squere);
    }

    Snap_rounding_rotation<base_rep> r;
    r(ms1,angle);
    r(ms2,angle);
    r(ms3,angle);
    r(ms4,angle);
    r(ms5,angle);
    r(ms6,angle);
    /*    _gt->rotate_point_2_object()(ms1,angle);
    _gt->rotate_point_2_object()(ms2,angle);
    _gt->rotate_point_2_object()(ms3,angle);
    _gt->rotate_point_2_object()(ms4,angle);
    _gt->rotate_point_2_object()(ms5,angle);
    _gt->rotate_point_2_object()(ms6,angle);*/

    // query
    Point_2 point_left,point_right,point_bot,point_top;

    point_left = small_x_point(ms1,ms2);
    point_left = small_x_point(point_left,ms3);
    point_left = small_x_point(point_left,ms4);
    point_left = small_x_point(point_left,ms5);
    point_left = small_x_point(point_left,ms6);

    point_right = big_x_point(ms1,ms2);
    point_right = big_x_point(point_right,ms3);
    point_right = big_x_point(point_right,ms4);
    point_right = big_x_point(point_right,ms5);
    point_right = big_x_point(point_right,ms6);

    point_bot = small_y_point(ms1,ms2);
    point_bot = small_y_point(point_bot,ms3);
    point_bot = small_y_point(point_bot,ms4);
    point_bot = small_y_point(point_bot,ms5);
    point_bot = small_y_point(point_bot,ms6);

    point_top = big_y_point(ms1,ms2);
    point_top = big_y_point(point_top,ms3);
    point_top = big_y_point(point_top,ms4);
    point_top = big_y_point(point_top,ms5);
    point_top = big_y_point(point_top,ms6);

    Iso_rectangle_2 rec(point_left,point_right,point_bot,point_top);

    return(rec);
  }

  friend class Snap_rounding_traits<base_rep>;
};

Bounding_box_of_minkowski_sum_2 bounding_box_of_minkowski_sum_2_object() const
    {return Bounding_box_of_minkowski_sum_2(this); }

 private:
  static std::map<const int,NT> angle_to_sines_appr;

  void init_angle_appr()
  {
    angle_to_sines_appr[0] = NT(0);
    angle_to_sines_appr[1] = NT(115) / NT(6613);
    angle_to_sines_appr[2] = NT(57) / NT(1625);
    angle_to_sines_appr[3] = NT(39) / NT(761);
    angle_to_sines_appr[4] = NT(29) / NT(421);
    angle_to_sines_appr[5] = NT(23) / NT(265);
    angle_to_sines_appr[6] = NT(19) / NT(181);
    angle_to_sines_appr[7] = NT(32) / NT(257);
    angle_to_sines_appr[8] = NT(129) / NT(929);
    angle_to_sines_appr[9] = NT(100) / NT(629);
    angle_to_sines_appr[10] = NT(92) / NT(533);
    angle_to_sines_appr[11] = NT(93) / NT(485);
    angle_to_sines_appr[12] = NT(76) / NT(365);
    angle_to_sines_appr[13] = NT(156) / NT(685);
    angle_to_sines_appr[14] = NT(205) / NT(853);
    angle_to_sines_appr[15] = NT(69) / NT(269);
    angle_to_sines_appr[16] = NT(7) / NT(25);
    angle_to_sines_appr[17] = NT(120) / NT(409);
    angle_to_sines_appr[18] = NT(57) / NT(185);
    angle_to_sines_appr[19] = NT(12) / NT(37);
    angle_to_sines_appr[20] = NT(51) / NT(149);
    angle_to_sines_appr[21] = NT(135) / NT(377);
    angle_to_sines_appr[22] = NT(372) / NT(997);
    angle_to_sines_appr[23] = NT(348) / NT(877);
    angle_to_sines_appr[24] = NT(231) / NT(569);
    angle_to_sines_appr[25] = NT(36) / NT(85);
    angle_to_sines_appr[26] = NT(39) / NT(89);
    angle_to_sines_appr[27] = NT(300) / NT(661);
    angle_to_sines_appr[28] = NT(8) / NT(17);
    angle_to_sines_appr[29] = NT(189) / NT(389);
    angle_to_sines_appr[30] = NT(451) / NT(901);
    angle_to_sines_appr[31] = NT(180) / NT(349);
    angle_to_sines_appr[32] = NT(28) / NT(53);
    angle_to_sines_appr[33] = NT(432) / NT(793);
    angle_to_sines_appr[34] = NT(161) / NT(289);
    angle_to_sines_appr[35] = NT(228) / NT(397);
    angle_to_sines_appr[36] = NT(504) / NT(865);
    angle_to_sines_appr[37] = NT(3) / NT(5);
    angle_to_sines_appr[38] = NT(580) / NT(941);
    angle_to_sines_appr[39] = NT(341) / NT(541);
    angle_to_sines_appr[40] = NT(88) / NT(137);
    angle_to_sines_appr[41] = NT(48) / NT(73);
    angle_to_sines_appr[42] = NT(65) / NT(97);
    angle_to_sines_appr[43] = NT(429) / NT(629);
    angle_to_sines_appr[44] = NT(555) / NT(797);
    angle_to_sines_appr[45] = NT(697) / NT(985);
    angle_to_sines_appr[46] = NT(572) / NT(797);
    angle_to_sines_appr[47] = NT(460) / NT(629);
    angle_to_sines_appr[48] = NT(72) / NT(97);
    angle_to_sines_appr[49] = NT(55) / NT(73);
    angle_to_sines_appr[50] = NT(105) / NT(137);
    angle_to_sines_appr[51] = NT(420) / NT(541);
    angle_to_sines_appr[52] = NT(741) / NT(941);
    angle_to_sines_appr[53] = NT(4) / NT(5);
    angle_to_sines_appr[54] = NT(703) / NT(865);
    angle_to_sines_appr[55] = NT(325) / NT(397);
    angle_to_sines_appr[56] = NT(240) / NT(289);
    angle_to_sines_appr[57] = NT(665) / NT(793);
    angle_to_sines_appr[58] = NT(45) / NT(53);
    angle_to_sines_appr[59] = NT(299) / NT(349);
    angle_to_sines_appr[60] = NT(780) / NT(901);
    angle_to_sines_appr[61] = NT(340) / NT(389);
    angle_to_sines_appr[62] = NT(15) / NT(17);
    angle_to_sines_appr[63] = NT(589) / NT(661);
    angle_to_sines_appr[64] = NT(80) / NT(89);
    angle_to_sines_appr[65] = NT(77) / NT(85);
    angle_to_sines_appr[66] = NT(520) / NT(569);
    angle_to_sines_appr[67] = NT(805) / NT(877);
    angle_to_sines_appr[68] = NT(925) / NT(997);
    angle_to_sines_appr[69] = NT(352) / NT(377);
    angle_to_sines_appr[70] = NT(140) / NT(149);
    angle_to_sines_appr[71] = NT(35) / NT(37);
    angle_to_sines_appr[72] = NT(176) / NT(185);
    angle_to_sines_appr[73] = NT(391) / NT(409);
    angle_to_sines_appr[74] = NT(24) / NT(25);
    angle_to_sines_appr[75] = NT(260) / NT(269);
    angle_to_sines_appr[76] = NT(828) / NT(853);
    angle_to_sines_appr[77] = NT(667) / NT(685);
    angle_to_sines_appr[78] = NT(357) / NT(365);
    angle_to_sines_appr[79] = NT(476) / NT(485);
    angle_to_sines_appr[80] = NT(525) / NT(533);
    angle_to_sines_appr[81] = NT(621) / NT(629);
    angle_to_sines_appr[82] = NT(920) / NT(929);
    angle_to_sines_appr[83] = NT(255) / NT(257);
    angle_to_sines_appr[84] = NT(180) / NT(181);
    angle_to_sines_appr[85] = NT(264) / NT(265);
    angle_to_sines_appr[86] = NT(420) / NT(421);
    angle_to_sines_appr[87] = NT(760) / NT(761);
    angle_to_sines_appr[88] = NT(1624) / NT(1625);
    angle_to_sines_appr[89] = NT(6612) / NT(6613);
    angle_to_sines_appr[90] = NT(1);
    }

};

template<class base_rep>
std::map<const int,typename base_rep::FT> Snap_rounding_traits<base_rep>::angle_to_sines_appr;

CGAL_END_NAMESPACE

#endif // CGAL_ISR_2_TRAITS_H
