#ifndef ARRANGEMENT_OF_SPHERES_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

//#include <CGAL/Root_of_2.h>
#include <CGAL/Root_of_traits.h>

#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_3/Filtered_sphere_line_intersection.h>

#include <CGAL/Tools/Label.h>


struct Sphere_key{
  static const int BL=-2, TR=-1, TEMP=-3;
  Sphere_key(): id_(-4){}
  Sphere_key(int i): id_(i){}
  int to_index() const {return id_;}
  bool is_valid() const {
    return id_ > -4;
  }
  bool operator==(const Sphere_key &o) const {
    return id_==o.id_;
  }
  bool operator!=(const Sphere_key &o) const {
    return id_!=o.id_;
  }
  bool operator<(const Sphere_key &o) const {
    return id_<o.id_;
  }
  bool operator>(const Sphere_key &o) const {
    return id_>o.id_;
  }
  bool is_input() const {
    return id_>=0;
  }
  unsigned int input_index() const {
    CGAL_precondition(is_input());
    return id_;
  }

  unsigned int internal_index() const {
    return id_+3;
  }

  std::ostream &write(std::ostream &out) const {
    if (id_==TEMP) out << "<T>";
    else if (id_ == TR) out << "<TR>";
    else if (id_ == BL) out << "<BL>";
    else out << "<" << id_ << ">";
    return out;
}

  static Sphere_key temp_key() {return Sphere_key(TEMP);}
  static Sphere_key bl_key() {return Sphere_key(BL);}
  static Sphere_key tr_key() {return Sphere_key(TR);}

  int id_;
};

inline std::ostream &operator<<(std::ostream &out, Sphere_key sk) {
  return sk.write(out);
}

struct Arrangement_of_spheres_traits_3 {
  typedef CGAL::Gmpq NT;
  typedef CGAL::Root_of_traits<NT>::RootOf_2 Quadratic_NT;
  typedef CGAL::Cartesian<NT> Geometric_kernel;
  typedef Geometric_kernel::Sphere_3 Sphere_3;
  typedef Geometric_kernel::Point_3 Point_3;
  typedef Geometric_kernel::Plane_3 Plane_3;
  typedef Geometric_kernel::Vector_3 Vector_3;
  typedef Geometric_kernel::Segment_3 Segment_3;
  typedef Geometric_kernel::Line_3 Line_3;
  typedef Geometric_kernel::Line_2 Line_2;
  typedef Geometric_kernel::Point_2 Point_2;
  typedef Geometric_kernel::Vector_2 Vector_2;
  typedef Geometric_kernel::Circle_2 Circle_2;
  typedef Geometric_kernel::Segment_2 Segment_2;
  

  typedef Sphere_key Key;
  


  /*typedef CGAL::Cartesian<double> DK;
    typedef DK::Point_2 DPoint;
    typedef DK::Sphere_3 DSphere;*/

  typedef Sphere_line_intersection<Geometric_kernel> Sphere_point_3;
  typedef Filtered_sphere_line_intersection<Geometric_kernel, 2> Event_point_3;


  struct Intersect_with_sweep {
    Intersect_with_sweep(NT z): z_(z){}
    typedef Circle_2 result_type;
    typedef Sphere_3 argument_type;

    Circle_2 operator()(Sphere_3 s) const{
      NT r2= s.squared_radius()-  CGAL::square(s.center().z()-z_);
      CGAL_assertion(r2>=0);
      Circle_2  c(Point_2(s.center().x(), s.center().y()), r2);
      return c;
    }

    NT z_;
  };

  Intersect_with_sweep intersect_with_sweep_object(NT z) const {
    return Intersect_with_sweep(z);
  }

  struct Intersects_with_sweep {
    Intersects_with_sweep(NT z): z_(z){}
    typedef bool result_type;
    typedef Sphere_3 argument_type;

    bool operator()(Sphere_3 s) const{
      NT r2= s.squared_radius()-  CGAL::square(s.center().z()-z_);
      return r2 >=0;
    }

    NT z_;
  };

  Intersects_with_sweep intersects_with_sweep_object(NT z) const {
    return Intersects_with_sweep(z);
  }

  /*struct Sphere_location{
    enum Location {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, IN_BIT=16, OUT_BIT=32 };
    Sphere_location(Event_point_3 sp): sp_(sp){
      CGAL_precondition(sp.line().to_vector().x()==0);
      CGAL_precondition(sp.line().to_vector().y()==0);
      CGAL_precondition(sp.line().to_vector().z()==1);
    }

    int operator()(Sphere_3 s) const {
      int r=0;
      Event_point_3 a(s, sp_.line());
      if (a.is_valid()) {
	Event_point_3 b(s, sp_.line().opposite());
	CGAL::Comparison_result ca= a.compare(sp_, 2);
	CGAL::Comparison_result cb= b.compare(sp_, 2);
	
	if (ca == CGAL::SMALLER && cb == CGAL::LARGER) {
	  r |= IN_BIT;
	} else if (ca == CGAL::EQUAL || cb == CGAL::EQUAL) {
	  r |= IN_BIT;
	  r |= OUT_BIT;
	} else {
	  r |= OUT_BIT;
	}
	//std::cout << ca << " " << cb << std::endl;
      } else {
	r |= OUT_BIT;
      }
      Plane_3 lrp(Point_3(s.center()), Vector_3(1, 0, 0));
      Plane_3 tbp(Point_3(s.center()), Vector_3(0, 1, 0));
      CGAL::Oriented_side xo= oriented_side(lrp, sp_);
      CGAL::Oriented_side yo= oriented_side(tbp, sp_);
      if (xo!= CGAL::ON_NEGATIVE_SIDE) {
	r |= R_BIT;
      } 
      if (xo != CGAL::ON_POSITIVE_SIDE){
	r |= L_BIT;
      }
      if (yo != CGAL::ON_NEGATIVE_SIDE) {
	r |= T_BIT;
      } 
      if (yo != CGAL::ON_POSITIVE_SIDE) {
	r |= B_BIT;
      }
      return r;
    }

    static std::string decode(int i) {
      std::string r;
      if (i & IN_BIT) r += "I";
      if (i & L_BIT) r += "L";
      if (i & R_BIT) r += "R";
      if (i & T_BIT) r += "T";
      if (i & B_BIT) r += "B";
      return r;
    }

    Event_point_3 sp_;
  };
  
  Sphere_location sphere_location_object(Event_point_3 sp) const {
    return Sphere_location(sp);
    }*/

  Geometric_kernel geometric_kernel_object() const {
    return Geometric_kernel();
  }

  /*template <class C>
  static CGAL::Bounded_side bounded_side(Sphere_3 s,
					 const Sphere_point_3 &p) {
    Sphere_point_3 ls(s, p.line(), true);
    if (!ls.is_finite()) 
      return CGAL::ON_UNBOUNDED_SIDE;
    Sphere_point_3 rs(s, p.line(), false);
    if (ls > rs) std::swap(ls,rs);
    if (ls < p && rs > p) return CGAL::ON_BOUNDED_SIDE;
    else if (ls > p || rs < p ) return CGAL::ON_UNBOUNDED_SIDE;
  else return CGAL::ON_BOUNDARY;
  }*/



  /*static inline Point_3 closest_point(Line_3 l, Point_3 p) {
    NT t= (l.to_vector()*(s.center()-CGAL::ORIGIN) 
	   - l.to_vector()*(l.point()-CGAL::ORIGIN))/(l.to_vector()*l.to_vector());
    Point_3 cp= l.point()+t*l.to_vector();
    return cp;
    }*/
};

#endif
