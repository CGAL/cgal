#ifndef SKIN_SURFACE_QUADRATIC_SURFACE_H
#define SKIN_SURFACE_QUADRATIC_SURFACE_H

#include <CGAL/Weighted_point.h>

CGAL_BEGIN_NAMESPACE

/////////////////////////////////////////////////
//                                             //
//  Virtual base class of a quadratic surface  //
//  This implementation is specific for        //
//  skin surfaces                              //
//                                             //
/////////////////////////////////////////////////
template < class K >
class Skin_surface_quadratic_surface_3 {
public:
  typedef typename K::Point_3             Point;
  typedef typename K::Vector_3            Vector;
  typedef typename K::Segment_3           Segment;
  typedef typename K::RT                  RT;
  typedef Weighted_point<Point, RT> Weighted_point;

  Skin_surface_quadratic_surface_3(RT eps=10e-10) : eps(eps) {
  }
  virtual ~Skin_surface_quadratic_surface_3() {};
  // Construct the intersection point with the segment (p0,p1)
  Point to_surface(Point const &p0, Point const &p1) {
    CGAL_assertion (Sign(value(p1) * value(p0)) != POSITIVE);
    RT sq_d = squared_distance(p0,p1);
    Point pp0=p0, pp1=p1, mid;

    if (value(p1) < value(p0)) {
      std::swap(pp0, pp1);
    }

    while (sq_d > eps) {
      mid = midpoint(pp0,pp1);
      if (value(mid) > 0) {
	pp1 = mid;
      } else {
	pp0 = mid;
      }
      sq_d /= 4;
    }
    return mid;
  };
  virtual Point to_surface(Point const &p, Vector const &v) = 0;
  inline Point to_surface(Segment const &s) {
    return to_surface(s.source(), s.target());
  }

  // Gradient descent
  virtual Point to_surface(Point const &p0) = 0;

  // compute the function value of p
  virtual RT value(Point const &p) const = 0;
  // Compute the gradient in p
  virtual Vector gradient(Point const &p) = 0; 
  // Compute the normal in p (normalized gradient)
  virtual Vector normal(Point const &p) = 0; 


  // return the dimension of the delaunay simplex:
  virtual int dimension() const = 0;

  // return a continuous density function on the skin surface:
  virtual RT sq_density(Point const &p) = 0;

  RT eps;
};

/////////////////////////////////////////////////
//                                             //
//  Sphere                                     //
//                                             //
//  orient*|x-p|^2/s - P == 0                  //
//  orient*|x-p|^2/s - P == 0                  //
//                                             //
//  p: center of the sphere                    //
//  P: squared radius of the sphere            //
//  orient: orientation of the sphere          //
//                                             //
/////////////////////////////////////////////////
template < class K >
class Skin_surface_sphere_3 : public Skin_surface_quadratic_surface_3<K> {
public:
  typedef Skin_surface_quadratic_surface_3<K> Parent; 
  typedef typename Parent::Point              Point;
  typedef typename Parent::Vector             Vector;
  typedef typename Parent::RT                 RT;
  typedef typename Parent::Weighted_point     Weighted_point;

  Skin_surface_sphere_3(Weighted_point wp, RT s, int orient, RT eps=10e-10)
    : Skin_surface_quadratic_surface_3<K>(eps),
      wp(wp), s(s), orient(orient) {
    assert((orient == -1) || (orient == 1));
  }
	
  RT value(Point const &x) const {
    return orient * (squared_distance(x, wp)/s - wp.weight());
  }

  Point to_surface(Point const &p, Vector const &v) {
    Vector pc = p-wp;
    RT sq_d = v*v;
    RT top = -pc*v/sq_d; 

    RT extr_val = squared_distance(p+top*v, wp) - s*wp.weight();
		
    if (extr_val > 0) {
      // 			std::cerr << "im. intersection[" << dimension() << "] "
      // 								<<  extr_val << "\n"; 
      return p;
    }
    RT d = sqrt(-extr_val/sq_d);

    // t should be in [0,1]
    RT t, t1;
    t = top + d; t1 = top - d;
    if (t1*t1 < t*t) t = t1;

    CGAL_assertion(abs(value(p + t*v)) < 0.001);
    return p + t*v;
  }
  Point to_surface(Point const &p) {
    Vector pc = p-wp;
    return wp + sqrt(s*wp.weight()/(pc*pc))*pc;
  }

  Vector normal(Point const &p) {
    Vector n = gradient(p);
    return n/sqrt(n*n);
  }
	
  Vector gradient(Point const &p) {
    return orient*(p-wp);
  }
	
  int dimension() const {
    if (orient == 1) return 0; else return 3;
  }

  RT sq_density(Point const &p) {
    assert(wp.weight() > 0);
    // #ifdef WRITE_DEBUG
    // 		if ( std::abs(s*wp.weight() - squared_distance(wp, p)) > 1e-10) {
    // 			std::cerr << __FILE__ << " l" << __LINE__ << " dist: "
    // 								<< std::abs(s*wp.weight() - squared_distance(wp, p))
    // 								<< std::endl;
    // 		}
    // #endif
    return s*wp.weight();
  }
private:
  Weighted_point wp;
  RT s;
  int orient;
};

/* *************************************************
 *                  Hyperboloid
 *
 * Point p.
 * Let y = (p-wp)*t/|t|
 *     x = |p-wp|^2 - y^2
 *
 * orient*((1-s) x^2 - s y^2 + s(1-s) W = 0
 * orient*((1-s)|p-wp|^2 - y^2 + s(1-s)W = 0
 *
 ***************************************************/

template < class K >
class Skin_surface_hyperboloid_3 : public Skin_surface_quadratic_surface_3<K> {
public:
  typedef Skin_surface_quadratic_surface_3<K> Parent; 
  typedef typename Parent::Point              Point;
  typedef typename Parent::Vector             Vector;
  typedef typename Parent::RT                 RT;
  typedef typename Parent::Weighted_point     Weighted_point;

  Skin_surface_hyperboloid_3(Weighted_point wp, Vector t, RT s, int orient,
    RT eps=10e-10)
    : Skin_surface_quadratic_surface_3<K>(eps), wp(wp), t(t), s(s), orient(orient) {
    assert((orient == -1) || (orient == 1));
    sq_t = t*t;
  }
	
  RT value(Point const &x) const {
    Vector dir = x-wp;
    RT tmp = dir*t;
    tmp = tmp*tmp/sq_t;

    return orient * (dir*dir/s - tmp/(s*(1-s))) + wp.weight();
  }
//   Point to_surface(Point const &p0, Point const &p1) {
//     assert(value(p0) * value(p1) <= 0);

//     Vector p0c = p0-wp;
//     Vector p1p0 = p1-p0;
//     RT sq_d = p1p0*p1p0;
		
//     RT p0_sym = p0c*t;
//     RT p1_sym = (p1-wp)*t;
//     RT d_sym = p0_sym-p1_sym;

//     RT den = ((1-s)*sq_d - d_sym*d_sym/sq_t);
//     RT top = -((1-s)*(p0c*p1p0) + p0_sym*d_sym/sq_t) / den;
//     RT extr_val =
//       orient*((1-s)*p0c*p0c-p0_sym*p0_sym/sq_t-s*(1-s)*wp.weight());
			
//     {
//       Point extr = p0+top*p1p0;
//       Vector dir = extr-wp;
//       RT tmp = dir*t; tmp *= tmp/sq_t;
//       extr_val = orient*(-tmp + (1-s)*dir*dir + s*(1-s)*wp.weight());
//     }
//     RT d = sqrt(-orient*extr_val/den);

//     RT t, t1;
//     t = top + d; t1 = top - d;
//     if ((2*t1-1)*(2*t1-1) < (2*t-1)*(2*t-1)) t = t1;

//     if (t < 0) {
//       std::cerr << "Hl[" << dimension() <<"] " <<  t << "\n";
//       return p0;
//     } else if (t > 1) {
//       std::cerr << "Hh[" << dimension() <<"] " <<  t << "\n";
//       return  p1;
//     } else {
//       assert (std::abs(value(p0 + t*p1p0)) < 0.001);
//       return p0 + t*p1p0;
//     }
//   }
  Point to_surface(Point const &p, Vector const &v){
    Vector pc = p-wp;
    Point p1 = p + v;
    RT sq_d = v*v;
		
    RT p_sym = pc*t;
    RT p1_sym = (p1-wp)*t;
    RT d_sym = p_sym-p1_sym;

    RT den = ((1-s)*sq_d - d_sym*d_sym/sq_t);
    RT top = -((1-s)*(pc*v) + p_sym*d_sym/sq_t) / den;
    RT extr_val =
      orient*((1-s)*pc*pc-p_sym*p_sym/sq_t-s*(1-s)*wp.weight());
			
    {
      Point extr = p+top*v;
      Vector dir = extr-wp;
      RT tmp = dir*t; tmp *= tmp/sq_t;
      extr_val = orient*(-tmp + (1-s)*dir*dir + s*(1-s)*wp.weight());
    }
    RT d = -orient*extr_val/den;
    if (d<0) return p;
    d = sqrt(d);

    RT t, t1;
    t = top + d; t1 = top - d;
    if (t1*t1 < t*t) t = t1;

    assert (std::abs(value(p + t*v)) < 0.001);
    return p + t*v;
  }
  Point to_surface(Point const &p0) {
    return to_surface(p0, normal(p0));
  }
  Vector normal(Point const &p) {
    Vector n = gradient(p);
    return orient*n/sqrt(n*n);
  }
  Vector gradient(Point const &p) {
    // -s x + (1-s) y 
    Vector v = p - wp;
    Vector vt = (v*t)/(t*t)*t;
    return orient*((1-s)*v - vt);
  }

  int dimension() const {
    if (orient == 1) return 1; else return 2;
  }

  RT sq_density(Point const &p) {
    Vector v = p - wp;
    RT vt = v*t;
    RT scale = 1-s;
    if (s < .5) {
      scale = (scale-s)/(scale*scale);
    } else {
      scale = (s-scale)/(s*s);
    }
    assert(scale >= 0);
		
    return squared_distance(wp, p) - scale*vt*vt/(t*t);
  }
private:
  Weighted_point wp;
  Vector t;
  RT s, sq_t;
  int orient;
	
};

CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_QUADRATIC_SURFACE_H
