#ifndef CEP_LEDA_RAT_GENERALIZED_PREDICATES_3_H
#define CEP_LEDA_RAT_GENERALIZED_PREDICATES_3_H

#include <CGAL/basic.h>

#include <CEP/Leda_rat_kernel/LEDA_RATKERNEL/d3_rat_support_functions.h>

// LEDA rational kernel generalized predicates ...

// 3d generalized predicates  ...

CGAL_BEGIN_NAMESPACE


class Predicate_leda_d3_rat_angle {
public:
  typedef Arity_tag< 3 > Arity;
  typedef CGAL::Angle    result_type;
  
  // we just return the sign of the dot product ...
  
  CGAL::Angle operator()(const leda_d3_rat_point& p1,const leda_d3_rat_point& p2,
                         const leda_d3_rat_point& p3)
  {
    // get vectors from p2-p1 and from p2-p3
    leda_rat_vector v1 = p1-p2;
    leda_rat_vector v2 = p3-p2;
  
    leda_rational s_prod = v1*v2;
  
    if (s_prod == 0) return CGAL::RIGHT;
    if (s_prod >  0) return CGAL::ACUTE;
    return CGAL::OBTUSE;
  }
};


class Predicate_leda_d3_rat_equal {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
    return ( p1 == p2 );
  }

  bool operator()(const leda_rat_vector& v1, const leda_rat_vector& v2) const
  {
    return ( v1 == v2 );
  }
  
  bool operator()(LEDA_NAMESPACE_NAME::rat_direction& d1, LEDA_NAMESPACE_NAME::rat_direction& d2) const
  {
    leda_rat_vector v1 = d1.get_vector();
    leda_rat_vector v2 = d2.get_vector();
    
    CGAL_precondition( (v1.dim()==3) && (v2.dim()==3));
    
    leda_rational xc1 = v1[0];
    leda_rational xc2 = v2[0];
    leda_rational yc1 = v1[1];
    leda_rational yc2 = v2[1];
    leda_rational zc1 = v1[2];
    leda_rational zc2 = v2[2];    
    
    leda_rational q1,q2,q3;
    
    // equal directions ???    
    if (xc2==0) {
      if (xc1!=0) return false;
      else q1 = 0;
    }
    else q1 = xc1/xc2;
    
    if (yc2==0) {
      if (yc1!=0) return false;
      else q2 = 0;
    }
    else q2 = yc1/yc2;
    
    if (zc2==0) {
      if (zc1!=0) return false;
      else q3 = 0;
    }
    else q3 = zc1/zc2;     
    
    if (q1==q2 && q1==q3) return true;
    return false;
  }
  
  bool operator()(const leda_d3_rat_line& l1, const leda_d3_rat_line& l2) const
  {
    if (! (l1 == l2)) return false;
    
    // same direction ???
    leda_rat_vector v1 = l1.to_vector();
    leda_rat_vector v2 = l2.to_vector();
    
    // compare directions ...
    LEDA_NAMESPACE_NAME::rat_direction d1(v1);
    LEDA_NAMESPACE_NAME::rat_direction d2(v2);
    return this->operator()(d1,d2);
  }    
  
  bool operator()(const leda_d3_rat_plane& p1, const leda_d3_rat_plane& p2) const
  {
    return (p1 == p2);
  }
  
  bool operator()(const leda_d3_rat_ray& r1, const leda_d3_rat_ray& r2) const
  {
    if (r1.point1() != r2.point1()) return false;
    
    // compare directions ...
    leda_rat_vector v1 = r1.to_vector();
    leda_rat_vector v2 = r2.to_vector();
    
    // compare directions ...
    LEDA_NAMESPACE_NAME::rat_direction d1(v1);
    LEDA_NAMESPACE_NAME::rat_direction d2(v2);
    return this->operator()(d1,d2);    
  }  

  bool operator()(const leda_d3_rat_segment& s1, const leda_d3_rat_segment& s2) const
  {
    return (s1 == s2);
  }  

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s1,
                  const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s2) const
  {
    return (s1 == s2);
  }
#else
  bool operator()(const leda_d3_rat_sphere& s1, const leda_d3_rat_sphere& s2) const
  {
    if (s1.is_degenerate() && s2.is_degenerate()) return true;
    leda_d3_rat_point a1 = s1.point1();
    leda_d3_rat_point a2 = s1.point2();
    leda_d3_rat_point a3 = s1.point3();
    leda_d3_rat_point a4 = s1.point4();
    
    if (s2.contains(a1) && s2.contains(a2) && s2.contains(a3) && s2.contains(a4)) return true;
    return false;
  }   
#endif

  bool operator()(const leda_d3_rat_triangle& t1, const leda_d3_rat_triangle& t2) const
  {
    return (t1 == t2);
  }  
  
  bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s1,
                  const LEDA_NAMESPACE_NAME::d3_rat_simplex& s2) const
  {
    return (s1 == s2);
  }   
  
  bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic1,
                  const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic2) const
  {
     leda_d3_rat_point lower_left1 = ic1.vertex(0);
     leda_d3_rat_point lower_left2 = ic2.vertex(0);    
     leda_d3_rat_point upper_right1 = ic1.vertex(7);
     leda_d3_rat_point upper_right2 = ic2.vertex(7);  
     
     if (lower_left1==lower_left2 && upper_right1==upper_right2) return true;
     return false;       
  }     
   
};

class Predicate_leda_d3_rat_equal_x {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return true;
     if (leda_d3_rat_point::cmp_x(p1,p2) == 0) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_equal_y {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return true;   
     if (leda_d3_rat_point::cmp_y(p1,p2) == 0) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_equal_z {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return true;   
     if (leda_d3_rat_point::cmp_z(p1,p2) == 0) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_equal_xy {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return true;   
     if (leda_d3_rat_point::cmp_x(p1,p2) == 0 && leda_d3_rat_point::cmp_y(p1,p2) == 0) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_equal_xyz {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (LEDA_NAMESPACE_NAME::identical(p1,p2)) return true;   
     if (leda_d3_rat_point::cmp_xyz(p1,p2) == 0) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_less_x {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (leda_d3_rat_point::cmp_x(p1,p2) == -1) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_less_y {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (leda_d3_rat_point::cmp_y(p1,p2) == -1) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_less_z {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (leda_d3_rat_point::cmp_z(p1,p2) == -1) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_less_xy {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     int cmx = leda_d3_rat_point::cmp_x(p1,p2);
     
     if (cmx == -1) return true;
     else {
      if (cmx==0 && leda_d3_rat_point::cmp_y(p1,p2) == -1) return true;
     }
     return false;
   }
};

class Predicate_leda_d3_rat_less_xyz {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
   {
     if (leda_d3_rat_point::cmp_xyz(p1,p2) == -1) return true;
     return false;   
   }
};

class Predicate_leda_d3_rat_compare_x {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
     return ( (Comparison_result) leda_d3_rat_point::cmp_x(p1,p2));     
  }
};

class Predicate_leda_d3_rat_compare_y {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
     return ( (Comparison_result) leda_d3_rat_point::cmp_y(p1,p2));     
  }
};

class Predicate_leda_d3_rat_compare_z {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
     return ( (Comparison_result) leda_d3_rat_point::cmp_z(p1,p2));     
  }
};

class Predicate_leda_d3_rat_compare_xy {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
     int cmx = leda_d3_rat_point::cmp_x(p1,p2);    
     if (cmx != 0) return  ((Comparison_result) cmx);
     int cmy = leda_d3_rat_point::cmp_y(p1,p2);
          
     return ( (Comparison_result) cmy);     
  }
};

class Predicate_leda_d3_rat_compare_xyz {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2) const
  {
     return ( (Comparison_result) leda_d3_rat_point::cmp_xyz(p1,p2));     
  }
};


// filter this later ...
class Predicate_leda_d3_rat_less_signed_distance_to_plane {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_plane& pl,
                  const leda_d3_rat_point& p1, 
	          const leda_d3_rat_point& p2) const
  {
     int ori1 = pl.side_of(p1);
     int ori2 = pl.side_of(p2);
  
     leda_rational d1 = (pl.sqr_dist(p1))*ori1;
     leda_rational d2 = (pl.sqr_dist(p2))*ori2;     
     return (d1 < d2);
  }   
};



class Predicate_leda_d3_rat_less_distance_to_point {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                   const leda_d3_rat_point& p3) const
   {
     int res = LEDA_NAMESPACE_NAME::cmp_distances(p1,p2,p1,p3);
     
     if (res==-1) return true;
     return false;
   }
};

class Predicate_leda_d3_rat_compare_distance {
public:
  typedef Arity_tag< 3 > Arity;
  typedef Comparison_result       result_type;

  Comparison_result operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                               const leda_d3_rat_point& p3) const
  {
     return ( (Comparison_result) LEDA_NAMESPACE_NAME::cmp_distances(p1,p2,p1,p3));     
  }  
};


class Predicate_leda_d3_rat_collinear {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                   const leda_d3_rat_point& p3) const
   {
      return LEDA_NAMESPACE_NAME::collinear(p1,p2,p3);
   }
};

class Predicate_leda_d3_rat_coplanar {
public:
  typedef Arity_tag< 4 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2,
                   const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
   {
      return LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p4);   
   }
};


class Predicate_leda_d3_rat_orientation {
public:
  typedef Arity_tag< 4 > Arity;
  typedef Orientation       result_type;

  Orientation operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                         const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
  {
     return (Orientation) LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,p4);
  }
};

class Predicate_leda_d3_rat_coplanar_orientation {
public:
  typedef Orientation       result_type;
  typedef Arity_tag< 3 >          Arity;

  Orientation operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                         const leda_d3_rat_point& p3) const
  {
     if (LEDA_NAMESPACE_NAME::collinear(p1,p2,p3)) return CGAL::COLLINEAR;
     
     // return 2d orientation ...
     // find projection plane :
     Orientation ori_xy = (Orientation) (LEDA_NAMESPACE_NAME::orientation_xy(p1,p2,p3));
     if ( ori_xy != CGAL::COLLINEAR ) return ori_xy;
     
     Orientation ori_xz = (Orientation) (LEDA_NAMESPACE_NAME::orientation_xz(p1,p2,p3));
     if ( ori_xz != CGAL::COLLINEAR ) return ori_xz;
     
     return (Orientation) (LEDA_NAMESPACE_NAME::orientation_yz(p1,p2,p3));       
  }

  Orientation operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                         const leda_d3_rat_point& p3, const leda_d3_rat_point& p4) const
  {
     // Preconditions:
     CGAL_precondition(! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3));
     CGAL_precondition(LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p4));
     
     // find a plane; project and compute orientation ...
     Orientation ori_xy = (Orientation) (LEDA_NAMESPACE_NAME::orientation_xy(p1,p2,p3));

    if ( ori_xy != CGAL::COLLINEAR )
      // the projection in (x,y) - plane is OK
      return Orientation( ori_xy * LEDA_NAMESPACE_NAME::orientation_xy(p1,p2,p4) );

    // we have to project onto another plane :

    if ((leda_d3_rat_point::cmp_x(p1,p2) != 0) || (leda_d3_rat_point::cmp_x(p1,p3) != 0))
    {
      // projection into (x,z)-plane is ok
      return Orientation ( LEDA_NAMESPACE_NAME::orientation_xz(p1,p2,p3) * 
                           LEDA_NAMESPACE_NAME::orientation_xz(p1,p2,p4) ); 
    }
    
    // projection into (y,z)-plane
    return Orientation ( LEDA_NAMESPACE_NAME::orientation_yz(p1,p2,p3) * 
                         LEDA_NAMESPACE_NAME::orientation_yz(p1,p2,p4) );
         
  }
};

class Predicate_leda_d3_rat_coplanar_side_of_bounded_circle {
public:
  typedef Arity_tag< 4 > Arity;
  typedef Bounded_side       result_type;

  Bounded_side operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                          const leda_d3_rat_point& p3, const leda_d3_rat_point& test) const
  {
    // Preconditions:
    CGAL_precondition(! LEDA_NAMESPACE_NAME::collinear(p1,p2,p3));
    CGAL_precondition(LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,test));

    leda_d3_rat_point t1(0,0,0), t2(1,0,0), t3(0,1,0), t4(0,0,1);     
    
    // find a 4. point for the sphere ...
    int ori;
    leda_d3_rat_point p4 = ((ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,t1)) != 0) ? t1:
              ((ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,t2)) != 0) ? t2:
              ((ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,t3)) != 0) ? t3: 
              ((ori = LEDA_NAMESPACE_NAME::orientation(p1,p2,p3,t4)) != 0) ? t4: t4;
    
    // now region of sphere test:
    int res = (ori * LEDA_NAMESPACE_NAME::side_of_sphere(p1,p2,p3,p4,test));  
    if (res==0) return CGAL::ON_BOUNDARY;
    if (res==-1) return CGAL::ON_UNBOUNDED_SIDE;
    return CGAL::ON_BOUNDED_SIDE;        
  }
};

class Predicate_leda_d3_rat_side_of_oriented_sphere {
public:
  typedef Arity_tag< 5 > Arity;
  typedef Oriented_side       result_type;

  Oriented_side operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                           const leda_d3_rat_point& p3, const leda_d3_rat_point& p4,
			   const leda_d3_rat_point& test) const
  {
     return (Oriented_side) LEDA_NAMESPACE_NAME::side_of_sphere(p1,p2,p3,p4,test);
  }
};


class Predicate_leda_d3_rat_side_of_bounded_sphere {
public:
  typedef Bounded_side       result_type;
  typedef Arity_tag< 5 >     Arity;

  Bounded_side operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                          const leda_d3_rat_point& p3, const leda_d3_rat_point& p4,
			  const leda_d3_rat_point& test) const
  {
     int res = LEDA_NAMESPACE_NAME::region_of_sphere(p1,p2,p3,p4,test);
     if (res==0) return CGAL::ON_BOUNDARY;
     if (res==-1) return CGAL::ON_UNBOUNDED_SIDE;
     return CGAL::ON_BOUNDED_SIDE;     
  }

  // center of sphere must be in plane through p1,p2,p3 ...
  Bounded_side operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                          const leda_d3_rat_point& p3, const leda_d3_rat_point& test) const
  {
     leda_d3_rat_point m = leda_support::construct_circle_center_3(p1,p2,p3);
     int res = LEDA_NAMESPACE_NAME::cmp_distances(m,p1,m,test);
     if (res==0) return CGAL::ON_BOUNDARY;
     if (res==-1) return CGAL::ON_UNBOUNDED_SIDE;
     return CGAL::ON_BOUNDED_SIDE;     
  }
  
  Bounded_side operator()(const leda_d3_rat_point& p1, const leda_d3_rat_point& p2, 
                          const leda_d3_rat_point& test) const
  {
     // compute midpoint of sphere with diamater p1p2; then compare distances
     // we could filter this later ...
     leda_d3_rat_point m = LEDA_NAMESPACE_NAME::midpoint(p1,p2);
     int res = LEDA_NAMESPACE_NAME::cmp_distances(m,p1,m,test);
     if (res==0) return CGAL::ON_BOUNDARY;
     if (res==-1) return CGAL::ON_UNBOUNDED_SIDE;
     return CGAL::ON_BOUNDED_SIDE;
  } 

};

class Predicate_leda_d3_rat_is_degenerate {
public:
  typedef Arity_tag< 1 > Arity;
  typedef bool       result_type;

  bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& c) const
  {
    // are the vertices collinear ???
    return LEDA_NAMESPACE_NAME::collinear(c.vertex(0),c.vertex(5),c.vertex(7));
  }
  
  bool operator()(const leda_d3_rat_line& l) const
  {
    return (l.point1() == l.point2());
  }
  
  bool operator()(const leda_d3_rat_plane& pl) const
  {
#if (__LEDA__ >= 440)      
    return ((pl.A() == 0) && (pl.B() == 0) && (pl.C() == 0));
#else
    leda_rat_vector n = pl.normal();
    return ((n.X() == 0) && (n.Y() == 0) && (n.Z() == 0));
#endif    
  }
  
  bool operator()(const leda_d3_rat_ray& r) const
  {
    return (r.point1() == r.point2());
  }    
  
  bool operator()(const leda_d3_rat_segment& seg) const
  {
    return (seg.source() == seg.target());
  }

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& sp) const
  {
    return (sp.squared_radius() == 0);
  }
#else  
  bool operator()(const leda_d3_rat_sphere& sp) const
  {
    return (sp.is_degenerate());
  }
#endif  
  
  bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s) const
  {
    return s.is_degenerate();
  }
  
  bool operator()(const leda_d3_rat_triangle& t) const
  {
    return t.is_degenerate();  
  }            
};

class Predicate_leda_d3_rat_has_on {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_line& l, const leda_d3_rat_point& p) const
   {
     return l.contains(p);
   }
   
   bool operator()(const leda_d3_rat_ray& r, const leda_d3_rat_point& p) const
   {
     return r.contains(p);
   }   
   
   bool operator()(const leda_d3_rat_segment& s, const leda_d3_rat_point& p) const
   {
     return s.contains(p);  
   }     
   
   bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
   {
     return pl.contains(p);  
   }     

   bool operator()(const leda_d3_rat_triangle& t, const leda_d3_rat_point& p) const
   {
     // coplanar ?
     leda_d3_rat_point p1 = t.point1();
     leda_d3_rat_point p2 = t.point2();
     leda_d3_rat_point p3 = t.point3();
      
     if (! LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p)) return false;
     
     // triangle and point are coplanar ...
     if (t.is_degenerate()){
        if (p1==p2 && p2==p3){ // point
	  if (p==p1) return true;
	  return false;
	}
	// segment
	leda_d3_rat_segment seg(p1,p2);
	if (! seg.contains(p3)){  // we need p3 ...
	  seg = leda_d3_rat_segment(p1,p3);
	  if (! seg.contains(p2)) seg = leda_d3_rat_segment(p2,p3);
	}
	
	return seg.contains(p);
     }
     // triangle is not degenerate; that means that one of the following projections
     // is a (non-degenerate) triangle 
     leda_rat_triangle tp = leda_support::project_xy(t); 
     leda_rat_point    pp = p.project_xy();
     if (! tp.is_degenerate()){
        LEDA_NAMESPACE_NAME::region_kind rk = tp.region_of(pp);
	if (rk == LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) return false;
	return true;
     }
     
     tp = leda_support::project_yz(t);
     pp = p.project_yz();
     
     if (! tp.is_degenerate()){
        LEDA_NAMESPACE_NAME::region_kind rk = tp.region_of(pp);
	if (rk == LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) return false;
	return true;
     }
     
     tp = leda_support::project_xz(t);
     pp = p.project_xz();
     LEDA_NAMESPACE_NAME::region_kind rk = tp.region_of(pp);
     if (rk == LEDA_NAMESPACE_NAME::UNBOUNDED_REGION) return false;
     return true;               
   }   
};

class Predicate_leda_d3_rat_has_on_bounded_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     if (d < sq) return true;
     return false;
      
  }
#else
   bool operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     return s.inside(p);
   }
#endif   

   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
     // inside returns true if p is inside or on s ...
     
     bool inside = s.in_simplex(p);
     if (! inside) return false;
     
     // is the point ON the simplex ???
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();   
     
     if (LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p) ||  LEDA_NAMESPACE_NAME::coplanar(p1,p2,p4,p) ||
         LEDA_NAMESPACE_NAME::coplanar(p2,p3,p4,p) ||  LEDA_NAMESPACE_NAME::coplanar(p3,p1,p4,p)) return false; 
     
     return true;
   }   
   
   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic, const leda_d3_rat_point& p) const
   {
     leda_d3_rat_point min = ic.vertex(0);
     leda_d3_rat_point max = ic.vertex(7);
     
     // compare with x/y/z coordinates ...
     if (leda_d3_rat_point::cmp_x(p,min) != 1 || leda_d3_rat_point::cmp_x(p,max) != -1) return false;
     if (leda_d3_rat_point::cmp_y(p,min) != 1 || leda_d3_rat_point::cmp_y(p,max) != -1) return false;
     if (leda_d3_rat_point::cmp_z(p,min) != 1 || leda_d3_rat_point::cmp_z(p,max) != -1) return false;     
          
     return true;     
   }      
};

class Predicate_leda_d3_rat_has_on_unbounded_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     if (d > sq) return true;
     return false;
      
  }
#else
   bool operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     return s.outside(p);
   }
#endif

   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
     bool inside = s.in_simplex(p);
     if (inside) return false;
     
     // is the point ON the simplex ???
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();   
     
     if (LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p) ||  LEDA_NAMESPACE_NAME::coplanar(p1,p2,p4,p) ||
         LEDA_NAMESPACE_NAME::coplanar(p2,p3,p4,p) ||  LEDA_NAMESPACE_NAME::coplanar(p3,p1,p4,p)) return false; 
     
     return true;
   }   
   
   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic, const leda_d3_rat_point& p) const
   {
     leda_d3_rat_point min = ic.vertex(0);
     leda_d3_rat_point max = ic.vertex(7);
     
     // compare with x/y/z coordinates ...
     if (leda_d3_rat_point::cmp_x(p,min) != 1 || leda_d3_rat_point::cmp_x(p,max) != -1) return true;
     if (leda_d3_rat_point::cmp_y(p,min) != 1 || leda_d3_rat_point::cmp_y(p,max) != -1) return true;
     if (leda_d3_rat_point::cmp_z(p,min) != 1 || leda_d3_rat_point::cmp_z(p,max) != -1) return true;     
          
     return false;     
   }      
};

class Predicate_leda_d3_rat_has_on_boundary {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;
  
  // undocumented
  bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
  {
      return pl.contains(p);
  }  

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     if (d == sq) return true;
     return false;
      
  }
#else
   bool operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     return s.contains(p);
   }
#endif

   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
     // s.inside returns true if p is in or on the simplex !!!!
     bool inside = s.in_simplex(p);
     if (! inside) return false;
     
     // is the point really ON the simplex ???
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();   
     
     if (LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p) ||  LEDA_NAMESPACE_NAME::coplanar(p1,p2,p4,p) ||
         LEDA_NAMESPACE_NAME::coplanar(p2,p3,p4,p) ||  LEDA_NAMESPACE_NAME::coplanar(p3,p1,p4,p)) return true; 
     
     return false;
   }   
   
   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic, const leda_d3_rat_point& p) const
   {
     leda_d3_rat_point min = ic.vertex(0);
     leda_d3_rat_point max = ic.vertex(7);
     
     // compare with x/y/z coordinates ...
     // we need at least one zero in the cmp results;
     int res1 = leda_d3_rat_point::cmp_x(p,min);
     int res2 = leda_d3_rat_point::cmp_x(p,max);
     
     if (res1 == -1 || res2 == 1) return false;
     
     int res3 = leda_d3_rat_point::cmp_y(p,min);
     int res4 = leda_d3_rat_point::cmp_y(p,max);     
     
     if (res3 == -1 || res4 == 1) return false;
     
     int res5 = leda_d3_rat_point::cmp_z(p,min);
     int res6 = leda_d3_rat_point::cmp_z(p,max);        
     
     if (res5 == -1 || res6 == 1) return false;     
          
     // we need at least one zero in the comparison results ...	  
     if (res1==0 || res2==0 || res3==0 || res4==0 || res5==0 || res6==0) return true;  	  
	  
     return false;     
   }      
};

class Predicate_leda_d3_rat_has_on_positive_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
   {
     int res = pl.side_of(p);
     if (res == 1) return true;
     return false;
   }

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     CGAL::Orientation ori = s.orientation();  
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     bool res;
     if (d < sq) res = true;
     else res = false;
     
     // testen ...
     if (ori == CLOCKWISE) res = !res;
     
     return res; 
  }
#else
   bool operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     int res = LEDA_NAMESPACE_NAME::side_of_sphere(s.point1(),s.point2(),s.point3(),s.point4() ,p);
     if (res == 1) return true;
     return false;
   }
#endif
   
//  welche Seite von Tetrahedron ist positive, welche ist die negative ? Dok. ist etwas unklar ...   

   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
       // orientation * bounded side
       int ori = LEDA_NAMESPACE_NAME::orientation(s.point1(), s.point2(), s.point3(), s.point4());
       
       if (ori != 0) {
       
         Bounded_side reg = leda_support::bounded_side(s,p);
	 
	 if (((CGAL::Oriented_side) (reg*ori)) == ON_POSITIVE_SIDE) return true;
	 else return false;
       }
       
       return false;      
   }
   
};

class Predicate_leda_d3_rat_has_on_negative_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

   bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
   {
     int res = pl.side_of(p);
     if (res == -1) return true;
     return false;
   }

#if defined(CGAL_COMPATIBLE_SPHERES)
  bool operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     CGAL::Orientation ori = s.orientation();  
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     bool res;
     if (d > sq) res = true;
     else res = false;
     
     // testen ...
     if (ori == CLOCKWISE) res = !res;
     
     return res; 
  }
#else
   bool operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     int res = LEDA_NAMESPACE_NAME::side_of_sphere(s.point1(),s.point2(),s.point3(),s.point4() ,p);
     if (res == -1) return true;
     return false;
   }
#endif
   
//  welche Seite von Tetrahedron ist positive, welche ist die negative ? Dok. ist etwas unklar ...   

   bool operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
       // orientation * bounded side
       int ori = LEDA_NAMESPACE_NAME::orientation(s.point1(), s.point2(), s.point3(), s.point4());
       
       if (ori != 0) {
       
         Bounded_side reg = leda_support::bounded_side(s,p);
	 
	 if (((CGAL::Oriented_side) (reg*ori)) == ON_NEGATIVE_SIDE) return true;
	 else return false;
       }
       
       return false;      
   }
   
};

class Predicate_leda_d3_rat_oriented_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Oriented_side       result_type;

   Oriented_side operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_point& p) const
   {
     return (Oriented_side) (pl.side_of(p));
   }

#if defined(CGAL_COMPATIBLE_SPHERES)
  CGAL::Oriented_side operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s,
                                 const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     CGAL::Orientation ori = s.orientation();  
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     if (d == sq) return CGAL::ON_ORIENTED_BOUNDARY;
     
     bool res;
     if (d < sq) res = true;
     else res = false;
     
     // testen ...
     if (ori == CLOCKWISE) res = !res;
     
     if (res) return CGAL::ON_POSITIVE_SIDE;
     return CGAL::ON_NEGATIVE_SIDE; 
  }
#else   
   Oriented_side operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& p) const
   {
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();     
     return (Oriented_side) (LEDA_NAMESPACE_NAME::side_of_sphere(p1,p2,p3,p4,p));
   } 
#endif

//  welche Seite von Tetrahedron ist positive, welche ist die negative ? Dok. ist etwas unklar ...  
//  precondition seems to be that the tetrahedron is not degenerate ...
 
   Oriented_side operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
       // orientation * bounded side
       int ori = LEDA_NAMESPACE_NAME::orientation(s.point1(), s.point2(), s.point3(), s.point4());
       
       if (ori != 0) {
       
         Bounded_side reg = leda_support::bounded_side(s,p);
	 
	 return (CGAL::Oriented_side) (reg*ori);
       }
       
       return ON_ORIENTED_BOUNDARY;     
   }
  
};

class Predicate_leda_d3_rat_bounded_side {
public:
  typedef Arity_tag< 2 > Arity;
  typedef Bounded_side   result_type;

#if defined(CGAL_COMPATIBLE_SPHERES)
  Bounded_side operator()(const LEDA_NAMESPACE_NAME::cgal_d3_rat_sphere& s, const leda_d3_rat_point& p) const
  {
     leda_d3_rat_point center = s.center();
     leda_rational sq      = s.squared_radius();
     leda_rational d       = center.sqr_dist(p);
     
     if (d==sq) return CGAL::ON_BOUNDARY;
     if (d <sq) return CGAL::ON_BOUNDED_SIDE;
     return CGAL::ON_UNBOUNDED_SIDE;
  }
#else
   Bounded_side operator()(const leda_d3_rat_sphere& s, const leda_d3_rat_point& test) const
   {
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();     
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();     
       
     int res = LEDA_NAMESPACE_NAME::region_of_sphere(p1,p2,p3,p4,test);
     if (res==0) return CGAL::ON_BOUNDARY;
     if (res==-1) return CGAL::ON_UNBOUNDED_SIDE;
     return CGAL::ON_BOUNDED_SIDE;     
   }
#endif   
  
   Bounded_side operator()(const LEDA_NAMESPACE_NAME::d3_rat_simplex& s, const leda_d3_rat_point& p) const
   {
     // inside returns true if p is inside or on s ...
     
     bool inside = s.in_simplex(p);
     if (! inside) return CGAL::ON_UNBOUNDED_SIDE;
     
     // is the point ON the simplex ???
     leda_d3_rat_point p1 = s.point1();
     leda_d3_rat_point p2 = s.point2();
     leda_d3_rat_point p3 = s.point3();
     leda_d3_rat_point p4 = s.point4();   
     
     if (LEDA_NAMESPACE_NAME::coplanar(p1,p2,p3,p) ||  LEDA_NAMESPACE_NAME::coplanar(p1,p2,p4,p) ||
         LEDA_NAMESPACE_NAME::coplanar(p2,p3,p4,p) ||  LEDA_NAMESPACE_NAME::coplanar(p3,p1,p4,p)) return CGAL::ON_BOUNDARY; 
     
     return CGAL::ON_BOUNDED_SIDE;
   }   
   
   Bounded_side operator()(const LEDA_NAMESPACE_NAME::d3_rat_iso_cuboid& ic, const leda_d3_rat_point& p) const
   {
     leda_d3_rat_point min = ic.vertex(0);
     leda_d3_rat_point max = ic.vertex(7);
     
     // compare with x/y/z coordinates ...
     // we need at least one zero in the cmp results;
     int res1 = leda_d3_rat_point::cmp_x(p,min);
     int res2 = leda_d3_rat_point::cmp_x(p,max);
     
     if (res1 == -1 || res2 == 1) return CGAL::ON_UNBOUNDED_SIDE;
     
     int res3 = leda_d3_rat_point::cmp_y(p,min);
     int res4 = leda_d3_rat_point::cmp_y(p,max);     
     
     if (res3 == -1 || res4 == 1) return CGAL::ON_UNBOUNDED_SIDE;
     
     int res5 = leda_d3_rat_point::cmp_z(p,min);
     int res6 = leda_d3_rat_point::cmp_z(p,max);        
     
     if (res5 == -1 || res6 == 1) return CGAL::ON_UNBOUNDED_SIDE;     
          
     // we need at least one zero in the comparison results ...	  
     if (res1==0 || res2==0 || res3==0 || res4==0 || res5==0 || res6==0) return CGAL::ON_BOUNDARY;  	  
	  
     return CGAL::ON_BOUNDED_SIDE;         
   }        
     
};

class Predicate_leda_d3_rat_are_ordered_along_line {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
  {
      return leda_d3_rat_segment(p,r).contains(q);  
  }
};

class Predicate_leda_d3_rat_are_strictly_ordered_along_line {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
  {
      return (leda_d3_rat_segment(p,r).contains(q) && ( q != p ) && ( q != r ));   
  }
};

// these two predicates could probably be optimized further ....
// (because the collinearity check in contains is not necessary -
// collinearity is a precondition) 

class Predicate_leda_d3_rat_collinear_are_ordered_along_line {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
  {
      return leda_d3_rat_segment(p,r).contains(q);  
  }
};

class Predicate_leda_d3_rat_collinear_are_strictly_ordered_along_line {
public:
  typedef Arity_tag< 3 > Arity;
  typedef bool       result_type;

  bool operator()(const leda_d3_rat_point& p, const leda_d3_rat_point& q, const leda_d3_rat_point& r) const
  {
      return (leda_d3_rat_segment(p,r).contains(q) && ( q != p ) && ( q != r ));   
  }
};

class Predicate_leda_d3_rat_do_intersect {
public:
  typedef Arity_tag< 2 > Arity;
  typedef bool       result_type;

// intersection tests with 3d plane ...

  bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_plane& pl2) const
  {
     if (pl.parallel(pl2)) return false;
     
     // non-parallel planes always intersect ...
     return true;
  }
  
  bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_line& l) const
  {
    leda_d3_rat_point p1 = l.point1();
    leda_d3_rat_point p2 = l.point2();
  
    int s1 = pl.side_of(p1);
    int s2 = pl.side_of(p2);
    
    if (s1==0 && s2==0) return true;
    if (s1 != s2) return true;  
    
    // the points p1 and p2 are on one side of the plane ...
    // so compare the distance of them; if it is not equal, 
    // we have an intersection ...
    // later we should have a filtered distance comparison predicate
    
    if (pl.sqr_dist(p1) == pl.sqr_dist(p2)) return false; // no intersection
    
    return true;
  }
  
  bool operator()(const leda_d3_rat_line& l, const leda_d3_rat_plane& pl) const
  {
     return this->operator()(pl, l);
  }  

  bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_ray& r) const
  {
    leda_d3_rat_point p1 = r.point1(); // source
    leda_d3_rat_point p2 = r.point2(); // target
  
    int s1 = pl.side_of(p1);
    int s2 = pl.side_of(p2);
    
    if (s1==0 && s2==0) return true;
    if (s1 != s2) return true;  
    
    // source and other point are on same side ...
    leda_rational r1 = pl.sqr_dist(p1);
    leda_rational r2 = pl.sqr_dist(p2); 
    
    if (r1 <= r2) return false; // the ray is going away from the plane
    return true;  
  }  

  bool operator()(const leda_d3_rat_ray& r, const leda_d3_rat_plane& pl) const
  {
     return this->operator()(pl, r);  
  }  

  bool operator()(const leda_d3_rat_plane& pl, const leda_d3_rat_segment& s) const
  {
    int s1 = pl.side_of(s.source());
    int s2 = pl.side_of(s.target());
    
    if (s1==0 && s2==0) return true;
    if (s1 != s2) return true;
    return false;
  }  
  
  bool operator()(const leda_d3_rat_segment& s, const leda_d3_rat_plane& pl) const
  {
     return this->operator()(pl, s);    
  }  
  
};


CGAL_END_NAMESPACE

#endif






