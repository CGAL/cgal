// LEDA 3d rat-kernel specializations ...
// see corresponding CGAL file for original version

#ifndef CONVEX_HULL_PROJECTIVE_YZ_TRAITS_LEDA_RAT_2_H
#define CONVEX_HULL_PROJECTIVE_YZ_TRAITS_LEDA_RAT_2_H

#include <CGAL/Convex_hull_projective_yz_traits_2.h>
#include <CGAL/leda_rational.h>
#include <LEDA/d3_rat_point.h>


namespace CGAL {

template <>
class Less_xy_plane_yz_2<leda_d3_rat_point>  
{
public:
   typedef bool         result_type;
   typedef Arity_tag<2> Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.ycoord(), p.zcoord(), q.ycoord(), q.zcoord()) == SMALLER;
   }
};

template <>
class Equal_xy_plane_yz_2<leda_d3_rat_point> 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<2> Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.ycoord(), p.zcoord(), q.ycoord(), q.zcoord()) == EQUAL;
   }
};

template <>
class Less_yx_plane_yz_2<leda_d3_rat_point>  
{
public:
   typedef bool         result_type;
   typedef Arity_tag<2> Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.zcoord(), p.ycoord(), q.zcoord(), q.ycoord()) == SMALLER;
   }
};

template <>
class Left_turn_plane_yz_2<leda_d3_rat_point>  
{
public:
   typedef bool         result_type;
   typedef Arity_tag<3> Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
   { 
    return orientationC2(p.ycoord(), p.zcoord(), q.ycoord(), q.zcoord(), r.ycoord(), r.zcoord()) == LEFT_TURN;
   }
};

template <>
class Left_of_line_plane_yz_2<leda_d3_rat_point> 
{
public:
   Left_of_line_plane_yz_2(const leda_d3_rat_point& a, const leda_d3_rat_point& b):
      p_a(a), p_b(b) 
   { }

   bool 
   operator()(const leda_d3_rat_point& c) const
   {
      return orientationC2(p_a.ycoord(), p_a.zcoord(), p_b.ycoord(), p_b.zcoord(), c.ycoord(), c.zcoord()) ==
             LEFT_TURN; 
   }
private:
   leda_d3_rat_point p_a;
   leda_d3_rat_point p_b;
};

template <>
class Less_dist_to_line_plane_yz_2<leda_d3_rat_point> 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<4> Arity;

   bool
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q,
              const leda_d3_rat_point& r, const leda_d3_rat_point& s) const
   {
      Comparison_result
         res = cmp_signed_dist_to_lineC2(p.ycoord(), p.zcoord(), q.ycoord(), q.zcoord(),
                                         r.ycoord(), r.zcoord(), s.ycoord(), s.zcoord());
      if ( res == LARGER )
         return false;
      else if ( res == SMALLER )
         return true;
      else
         return compare_lexicographically_xyC2(r.ycoord(), r.zcoord(), s.ycoord(), s.zcoord())
             == SMALLER;
   }
};

template <>
class Less_rotate_ccw_plane_yz_2<leda_d3_rat_point> 
{
public:
   typedef bool         result_type;
   typedef Arity_tag<3> Arity;

   bool
   operator()(const leda_d3_rat_point& r, const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   {
      Orientation orient =
               orientationC2(r.ycoord(), r.zcoord(), p.ycoord(), p.zcoord(), q.ycoord(), q.zcoord());
      if ( orient ==  LEFT_TURN )
         return true;
      else if ( orient == RIGHT_TURN )
         return false;
      else
      {
         if (p.ycoord() == r.ycoord() && p.zcoord() == r.zcoord()) return false;
         if (q.ycoord() == r.ycoord() && q.zcoord() == r.zcoord()) return true;
         if (p.ycoord() == q.ycoord() && p.zcoord() == q.zcoord()) return false;
         return
            collinear_are_ordered_along_lineC2(r.ycoord(), r.zcoord(),
                                               q.ycoord(), q.zcoord(), p.ycoord(), p.zcoord());
      }
   }
};



template <>
class Convex_hull_projective_yz_traits_2<leda_d3_rat_point>  
{
public:
    typedef leda_d3_rat_point                  Point_3;
    typedef Point_3                            Point_2;
    typedef Less_xy_plane_yz_2<Point_3>        Less_xy_2;
    typedef Equal_xy_plane_yz_2<Point_3>       Equal_2;
    typedef Less_yx_plane_yz_2<Point_3>        Less_yx_2;
    typedef Left_turn_plane_yz_2<Point_3>      Leftturn_2;
    typedef Left_turn_plane_yz_2<Point_3>      Left_turn_2;    
    
    typedef Less_dist_to_line_plane_yz_2<Point_3>
                                                Less_signed_distance_to_line_2;
    typedef Less_rotate_ccw_plane_yz_2<Point_3> Less_rotate_ccw_plane_2;

    Less_xy_2
    less_xy_2_object() const
    {  return Less_xy_2(); }
        
    Equal_2
    equal_2_object() const
    {  return Equal_2(); }        

    Less_yx_2
    less_yx_2_object() const
    {  return Less_yx_2(); }

    Leftturn_2
    leftturn_2_object() const
    {  return Leftturn_2(); }
    
    Left_turn_2
    left_turn_2_object() const
    {  return Left_turn_2(); }    

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {  return Less_signed_distance_to_line_2(); }

    Less_rotate_ccw_plane_2
    less_rotate_ccw_plane_2_object() const
    {  return Less_rotate_ccw_plane_2(); }
};

}
#endif 
