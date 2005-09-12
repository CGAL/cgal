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

  // Construct the intersection point with the segment (p0,p1)
  virtual Point to_surface(Point const &p0, Point const &p1) = 0;
  virtual Point to_surface(Point const &p, Vector const &v) = 0;
  inline Point to_surface(Segment const &s) {
    return to_surface(s.source(), s.target());
  }
  virtual Point to_surface(Point const &p0) = 0;

  // Compute the normal in p:
  virtual Vector normal(Point const &p) = 0; 

  // compute the function value of p
  virtual RT value(Point const &p) = 0;

  // return the dimension of the delaunay simplex:
  virtual int dimension() = 0;

  // return a continuous density function on the skin surface:
  virtual RT sq_density(Point const &p) = 0;

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
  typedef typename K::Point_3             Point;
  typedef typename K::Vector_3            Vector;
  typedef typename K::RT                  RT;
  typedef Weighted_point<Point, RT> Weighted_point;

  Skin_surface_sphere_3(Weighted_point wp, RT s, int orient)
    : Skin_surface_quadratic_surface_3<K>(), wp(wp), s(s), orient(orient) {
    assert((orient == -1) || (orient == 1));
  }
	
  RT value(Point const &x) {
    return orient*(1-s)*(CGAL::squared_distance(x, wp) - s*wp.weight());
  }
  Point to_surface(Point const &p0, Point const &p1) {
    Vector p0c = p0-wp;
    Vector p0p1 = p1-p0;
    RT sq_d = p0p1*p0p1;
    RT top = -p0c*p0p1/sq_d; 

    RT extr_val = CGAL::squared_distance(p0+top*p0p1, wp) - s*wp.weight();
		
    if (extr_val > 0) {
      std::cerr << "im. intersection[" << dimension() << "] "
		<<  extr_val << "\n"; 
      return p0 + top*p0p1;
    }
    RT d = sqrt(-extr_val/sq_d);

    // t should be in [0,1]
    RT t, t1;
    t = top + d; t1 = top - d;
    if ((2*t1-1)*(2*t1-1) < (2*t-1)*(2*t-1)) t = t1;

    if (t < 0) {
      std::cerr << "Sl[" << dimension() <<"] " <<  t << "\n";
      return  p0;
    } else if (t > 1) {
      std::cerr << "Sh[" << dimension() <<"] " <<  t << "\n";
      return  p1;
    } else {
      // 			if (std::abs(value(p0 + t*p0p1)) >= 0.001) {
      // 				std::cerr << "VAL: " << value(p0 + t*p0p1) << std::endl;
      // 			}
      // 			assert (std::abs(value(p0 + t*p0p1)) < 0.001);
      return p0 + t*p0p1;
    }
  }
  Point to_surface(Point const &p, Vector const &v) {
    Vector pc = p-wp;
    RT sq_d = v*v;
    RT top = -pc*v/sq_d; 

    RT extr_val = CGAL::squared_distance(p+top*v, wp) - s*wp.weight();
		
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
    Vector n = orient*(p-wp);
    return n/sqrt(n*n);
  }
	
  int dimension() {
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
  typedef typename K::Point_3             Point;
  typedef typename K::Vector_3            Vector;
  typedef typename K::RT                  RT;
  typedef CGAL::Weighted_point<Point, RT> Weighted_point;

  Skin_surface_hyperboloid_3(Weighted_point wp, Vector t, RT s, int orient)
    : Skin_surface_quadratic_surface_3<K>(), wp(wp), t(t), s(s), orient(orient) {
    assert((orient == -1) || (orient == 1));
    if (orient == -1) this->wp = Weighted_point(wp.point(), -wp.weight());
    sq_t = t*t;
  }
	
  RT value(Point const &x) {
    Vector dir = x-wp;
    RT tmp = dir*t;
    tmp = tmp*tmp/sq_t;

    return orient * ((1-s)*dir*dir - tmp + s*(1-s)*wp.weight());
  }
  Point to_surface(Point const &p0, Point const &p1) {
    assert(value(p0) * value(p1) <= 0);

    Vector p0c = p0-wp;
    Vector p1p0 = p1-p0;
    RT sq_d = p1p0*p1p0;
		
    RT p0_sym = p0c*t;
    RT p1_sym = (p1-wp)*t;
    RT d_sym = p0_sym-p1_sym;

    RT den = ((1-s)*sq_d - d_sym*d_sym/sq_t);
    RT top = -((1-s)*(p0c*p1p0) + p0_sym*d_sym/sq_t) / den;
    RT extr_val =
      orient*((1-s)*p0c*p0c-p0_sym*p0_sym/sq_t-s*(1-s)*wp.weight());
			
    {
      Point extr = p0+top*p1p0;
      Vector dir = extr-wp;
      RT tmp = dir*t; tmp *= tmp/sq_t;
      extr_val = orient*(-tmp + (1-s)*dir*dir + s*(1-s)*wp.weight());
    }
    RT d = sqrt(-orient*extr_val/den);

    RT t, t1;
    t = top + d; t1 = top - d;
    if ((2*t1-1)*(2*t1-1) < (2*t-1)*(2*t-1)) t = t1;

    if (t < 0) {
      std::cerr << "Hl[" << dimension() <<"] " <<  t << "\n";
      return p0;
    } else if (t > 1) {
      std::cerr << "Hh[" << dimension() <<"] " <<  t << "\n";
      return  p1;
    } else {
      assert (std::abs(value(p0 + t*p1p0)) < 0.001);
      return p0 + t*p1p0;
    }
  }
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
    // 		Vector n = normal(p0);
    // 		RT val1, val0 = value(p0);
    // 		Point p1 = p0 - val0*n;
    // 		val1 = value(p1);
    // 		while (std::abs(val1) > 1e-10)
    // 		{
    // 			val0 = val1;
    // 			p1 = p1 - val0*n;
    // 			val1 = value(p1);
    // 		}
    // 		return p1;
  }
  Vector normal(Point const &p) {
    // -s x + (1-s) y 
    Vector v = p - wp;
    Vector vt = (v*t)/(t*t)*t;
    Vector n = (1-s)*v - vt;
    return orient*n/sqrt(n*n);
  }

  int dimension() {
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
		
    return CGAL::squared_distance(wp, p) - scale*vt*vt/(t*t);
  }
private:
  Weighted_point wp;
  Vector t;
  RT s, sq_t;
  int orient;
	
};

CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_QUADRATIC_SURFACE_H
