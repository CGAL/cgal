// LEDA 3d rat-kernel specializations ...
// see corresponding CGAL file for original version

#ifndef CONVEX_HULL_PROJECTIVE_XZ_TRAITS_LEDA_RAT_2_H
#define CONVEX_HULL_PROJECTIVE_XZ_TRAITS_LEDA_RAT_2_H


#include <CGAL/Convex_hull_projective_xz_traits_2.h>
#include <CGAL/leda_rational.h>
#include <LEDA/d3_rat_point.h>

namespace CGAL {

template <>
class Less_xy_plane_xz_2<leda_d3_rat_point>  
{
public:
   typedef bool           result_type;
   typedef Arity_tag<2>   Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.xcoord(), p.zcoord(), q.xcoord(), q.zcoord()) == SMALLER;
   }
};

template <>
class Equal_xy_plane_xz_2<leda_d3_rat_point>
{
public:
   typedef bool           result_type;
   typedef Arity_tag<2>   Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.xcoord(), p.zcoord(), q.xcoord(), q.zcoord()) == EQUAL;
   }
};

template <>
class Less_yx_plane_xz_2<leda_d3_rat_point>  
{
public:
   typedef bool           result_type;
   typedef Arity_tag<2>   Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   { 
      return 
        compare_lexicographically_xyC2(p.zcoord(), p.xcoord(), q.zcoord(), q.xcoord()) == SMALLER;
   }
};

template <>
class Left_turn_plane_xz_2<leda_d3_rat_point>  
{
public:
   typedef bool           result_type;
   typedef Arity_tag<2>   Arity;

   bool 
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
   { 
    return orientationC2(p.xcoord(), p.zcoord(), q.xcoord(), q.zcoord(), r.xcoord(), r.zcoord()) == LEFTTURN;
   }
};


template <>
class Less_dist_to_line_plane_xz_2<leda_d3_rat_point> 
{
public:
   typedef bool           result_type;
   typedef Arity_tag<4>   Arity;

   bool
   operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q,
              const leda_d3_rat_point& r, const leda_d3_rat_point& s) const
   {
      Comparison_result
         res = cmp_signed_dist_to_lineC2(p.xcoord(), p.zcoord(), q.xcoord(), q.zcoord(),
                                         r.xcoord(), r.zcoord(), s.xcoord(), s.zcoord());
      if ( res == LARGER )
         return false;
      else if ( res == SMALLER )
         return true;
      else
         return compare_lexicographically_xyC2(r.xcoord(), r.zcoord(), s.xcoord(), s.zcoord())
             == SMALLER;
   }
};


template <>
class Less_rotate_ccw_plane_xz_2<leda_d3_rat_point> 
{
public:

   typedef bool         result_type;
   typedef Arity_tag<3> Arity;

   bool
   operator()(const leda_d3_rat_point& r, const leda_d3_rat_point& p, const leda_d3_rat_point& q) const
   {
      Orientation orient =
               orientationC2(r.xcoord(), r.zcoord(), p.xcoord(), p.zcoord(), q.xcoord(), q.zcoord());
      if ( orient ==  LEFTTURN )
         return true;
      else if ( orient == RIGHTTURN )
         return false;
      else
      {
         if (p.xcoord() == r.xcoord() && p.zcoord() == r.zcoord()) return false;
         if (q.xcoord() == r.xcoord() && q.zcoord() == r.zcoord()) return true;
         if (p.xcoord() == q.xcoord() && p.zcoord() == q.zcoord()) return false;
         return
            collinear_are_ordered_along_lineC2(r.xcoord(), r.zcoord(),
                                               q.xcoord(), q.zcoord(), p.xcoord(), p.zcoord());
      }
   }

};



template <>
class Convex_hull_projective_xz_traits_2<leda_d3_rat_point>  
{
public:
    typedef leda_d3_rat_point                   Point_3;
    typedef Point_3                             Point_2;
    typedef Less_xy_plane_xz_2<Point_3>         Less_xy_2;
    typedef Equal_xy_plane_xz_2<Point_3>        Equal_2;        
    typedef Less_yx_plane_xz_2<Point_3>         Less_yx_2;
    typedef Left_turn_plane_xz_2<Point_3>       Leftturn_2;
    typedef Left_turn_plane_xz_2<Point_3>       Left_turn_2;    
    
    typedef Less_rotate_ccw_plane_xz_2<Point_3> Less_rotate_ccw_2;
    typedef Less_dist_to_line_plane_xz_2<Point_3> 
                                                Less_signed_distance_to_line_2;

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

    Less_rotate_ccw_2
    less_rotate_ccw_2_object() const
    {  return Less_rotate_ccw_2(); }

    Less_signed_distance_to_line_2
    less_signed_distance_to_line_2_object() const
    {  return Less_signed_distance_to_line_2(); }

};

}
#endif 

