#ifndef SPHERE_LINE_INTERSECTION_H
#define SPHERE_LINE_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
#include <CGAL/tags.h>
#include <CGAL/Point_3.h>
#include <CGAL/Sqrt_extension.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

//#include <CGAL/Arrangement_of_spheres_traits_3.h>

//using doubles for ease of implementation, may not be the fastest
/*
  cases to handle:
  - tangent line -- fall back to algebraic since interval is small
  - tiny sphere--fall back to algebraic since interval is small
  - missing line
  - vertical line-- asserted
   
  I want to make sure that the interval is open. Now it is not. 
*/

template <class K>
struct Sphere_line_intersection {
private:
  /*static typename K::Point_3 unproject(typename K::Point_2 p) {
    return typename K::Point_3(p.x(), p.y(),0);
  }
  static typename K::Vector_3 unproject(typename K::Vector_2 p) {
    return typename K::Vector_3(p.x(), p.y(),0);
  }
  static typename K::Line_3 unproject(typename K::Line_2 p) {
    return typename K::Line_3(unproject(p.point()), unproject(p.to_vector()));
    }*/
public:
  typedef K T;
  typedef typename K::FT NT;
  CGAL_AOS3_DEFINE_QUADRATIC;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Line_3 Line_3;
  typedef typename K::Sphere_3 Sphere_3;
  typedef CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;

  //enum Type {FIRST, SECOND, ONE, INVALID};
  
  Sphere_line_intersection(): l_(Point_3(0,0,0), Vector_3(0,0,0)),
			      has_exact_(false){
    CGAL_assertion(l_.is_degenerate());
    if (0) {
      double d= approximate_coordinate(0);
      d+= interval_coordinate(0).first;
      std::cout << d << std::endl;
    }
  }

  /*Sphere_line_intersection(NT n):  type_(false){
    initialize_const(n);
    }*/
  Sphere_line_intersection(typename T::Point_3 p3,
			   typename T::Line_3 l): s_(p3,0), 
						  l_(l), 
						  has_exact_(false){
    CGAL_exactness_precondition(l.has_on(p3));
    CGAL_postcondition(is_valid());
  }


  Sphere_line_intersection(typename T::Point_3 p3 ): s_(p3,0),
						     l_(p3, sweep_vector<Vector_3>()),
    has_exact_(false){
    CGAL_postcondition(is_valid());
   
  }
  // assume the line hits
  /*Sphere_line_intersection( typename T::Sphere_3 s, 
			    typename T::Line_3 l, CGAL::Tag_true ): s_(s), 
								    l_(l),
								    has_exact_(false){
								    }*/
 
  Sphere_line_intersection(typename T::Sphere_3 s, typename T::Line_3 l): s_(s), 
									  l_(l), 
									  has_exact_(false) {
    typename T::Point_3 cp= l_.projection(s_.center()); //typename T::closest_point(l, s.center());
    CGAL::Bounded_side bs= s_.bounded_side(cp);
    if (bs == CGAL::ON_UNBOUNDED_SIDE){
      s_= Sphere_3();
      l_= Line_3(Point_3(0,0,0), Vector_3(0,0,0));
      CGAL_assertion(l_.is_degenerate());
      CGAL_assertion(!is_valid());
    } else if (bs== CGAL::ON_BOUNDARY) {
      s_= Sphere_3(cp, 0);
    }
  }
  
  // -b +- sqrt(b^2-4ac) / 2a
  // a0= (-b/2a)
  // a1= +-1/2a
  // sqrt= b^2-4ac

  
  bool has_exact() const {
    return has_exact_;
  }

  Bbox_3 bbox() const {
    return Bbox_3(interval_coordinate(Coordinate_index::X()).first,
		  interval_coordinate(Coordinate_index::Y()).first,
		  interval_coordinate(Coordinate_index::Z()).first,
		  interval_coordinate(Coordinate_index::X()).second,
		  interval_coordinate(Coordinate_index::Y()).second,
		  interval_coordinate(Coordinate_index::Z()).second);
  }

  typename CGAL::Root_of_traits<NT>::RootOf_2 
  root_of_2(Coordinate_index ci) const {
    Vector_3 lp= l_.point()-CGAL::ORIGIN;
    Vector_3 vc= s_.center()-CGAL::ORIGIN;
    Vector_3 lv= l_.to_vector(); 
    NT a=lv*lv;
    NT b=2*lv*(lp-vc); //-2*lv* vc + 2*lv*lp;
    NT c=lp*lp + vc*vc-s_.squared_radius()-2*lp*vc;
    
    NT disc= b*b-4*a*c;
    CGAL_assertion(disc >= 0);
   
    CGAL_precondition(a!=0);
    
    int coord =sweep_coordinate().index();
    CGAL::Sign sn= CGAL::sign(lv[coord]);
    bool first;
    if (sn != CGAL::ZERO) {
      first= sn == CGAL::POSITIVE;
    } else {
      coord =plane_coordinate(0).index();
      sn = CGAL::sign(lv[coord]);
      if (sn != CGAL::ZERO) {
	first= sn == CGAL::POSITIVE;
      } else {
	coord =plane_coordinate(1).index();
	sn = CGAL::sign(lv[coord]);
	first = sn == CGAL::POSITIVE;
      }
    }
    bool rt= first && CGAL::sign(lv[coord]) == CGAL::POSITIVE 
      || !first && CGAL::sign(lv[coord]) == CGAL::NEGATIVE;
    CGAL_assertion(first && CGAL::sign(lv[coord]) == CGAL::POSITIVE
		   || !first && CGAL::sign(lv[coord]) == CGAL::NEGATIVE);
    /*|| !first && CGAL::sign(lv[C]) == CGAL::NEGATIVE;*/
    typename CGAL::Root_of_traits<NT>::RootOf_2 t=make_root_of_2(a,b,c,rt);
    
	
	/*if (first) {
	  CGAL_assertion(t*lv[C] <= lv[C] *CGAL::make_root_of_2(a,b,c, !rt));
	  } else {
	  CGAL_assertion(t*lv[C] >= lv[C] *CGAL::make_root_of_2(a,b,c, !rt)); 
	  }*/
    return lp[ci.index()] + lv[ci.index()] *t;
  }
  

  void set_has_exact(bool tf) const {
    if (tf && !has_exact_) {
      Vector_3 lp= l_.point()-CGAL::ORIGIN;
      Vector_3 vc= s_.center()-CGAL::ORIGIN;
      Vector_3 lv= l_.to_vector(); 
      NT a=lv*lv;
      NT b=2*lv*(lp-vc); //-2*lv* vc + 2*lv*lp;
      NT c=lp*lp + vc*vc-s_.squared_radius()-2*lp*vc;
      
      NT disc= b*b-4*a*c;
      CGAL_assertion(disc >= 0);
      if (a==0) {
	CGAL_LOG(Log::SOME,l_ << std::endl << s_ << std::endl);
      }
      /*if (a==0) { 
	CGAL_assertion(s_.squared_radius() ==0);
	exact_[0]= s_.center()[0];
	exact_[1]= s_.center()[1];
	exact_[2]= s_.center()[2];
	} else {*/
	CGAL_precondition(a!=0);
	
	int coord =sweep_coordinate().index();
	CGAL::Sign sn= CGAL::sign(lv[coord]);
	bool first;
	if (sn != CGAL::ZERO) {
	  first= sn == CGAL::POSITIVE;
	} else {
	  coord =plane_coordinate(0).index();
	  sn = CGAL::sign(lv[coord]);
	  if (sn != CGAL::ZERO) {
	    first= sn == CGAL::POSITIVE;
	  } else {
	    coord =plane_coordinate(1).index();
	    sn = CGAL::sign(lv[coord]);
	    first = sn == CGAL::POSITIVE;
	  }
	}
	bool rt= first && CGAL::sign(lv[coord]) == CGAL::POSITIVE 
	  || !first && CGAL::sign(lv[coord]) == CGAL::NEGATIVE;
	CGAL_assertion(first && CGAL::sign(lv[coord]) == CGAL::POSITIVE
		       || !first && CGAL::sign(lv[coord]) == CGAL::NEGATIVE);
	/*|| !first && CGAL::sign(lv[C]) == CGAL::NEGATIVE;*/
	Quadratic_NT t=CGAL_AOS3_CONSTRUCT_QUADRATIC(a,b,c,rt);
	
	
	/*if (first) {
	  CGAL_assertion(t*lv[C] <= lv[C] *CGAL::make_root_of_2(a,b,c, !rt));
	  } else {
	  CGAL_assertion(t*lv[C] >= lv[C] *CGAL::make_root_of_2(a,b,c, !rt)); 
	  }*/
	for (int i=0; i< 3; ++i ) {
	  exact_[i]= lp[i] + lv[i]*t;
	}
	has_exact_=true;
      }
    //}
  }
  
  CGAL_IS(valid, 
	  return !l_.is_degenerate());

  CGAL_GET(Sphere_3, sphere, return s_);
  CGAL_GET(Line_3, line, return l_);
 
  Sphere_line_intersection flip_on_sweep() const {
    NT sp[3]= {s_.center().x(), s_.center().y(), s_.center().z()};
    sp[sweep_coordinate().index()]= - sp[sweep_coordinate().index()];
    Sphere_3 s(Point_3(sp[0], sp[1], sp[2]), s_.squared_radius());
    NT lp[3]= {l_.point().x(), l_.point().y(), l_.point().z()};
    lp[sweep_coordinate().index()]= - lp[sweep_coordinate().index()];
    NT lv[3]= {l_.to_vector().x(), l_.to_vector().y(), l_.to_vector().z()};
    lv[sweep_coordinate().index()]= - lv[sweep_coordinate().index()];
    Line_3 l(Point_3(lp[0], lp[1], lp[2]),
	     Vector_3(lv[0], lv[1], lv[2]));
    Sphere_line_intersection ret(s,l);
    if (has_exact_) {
      ret.has_exact_=true;
      for (unsigned int i=0; i< 3; ++i) {
	ret.exact_[i]= -exact_[i];
      }
    }
    return ret;
  }

  CGAL::Comparison_result compare_on_line(const Point_3 &pt) const {
    CGAL_exactness_precondition(line().has_on(pt));
    for (unsigned int i=0; i< 3; ++i){
      CGAL::Comparison_result c= compare(pt,Coordinate_index(i));
      if (c != CGAL::EQUAL) {
	if (line().to_vector()[i] >0) return c;
	else return CGAL::Comparison_result(-c);
      }
    }
    return CGAL::EQUAL;
  }



  std::ostream &write(std::ostream &out) const {
    //out << "(" << s_ << ", " << l_ << ", (" << lb_ << "..." << ub_ << "))";
    if (!is_valid()) {
      out << "INVALID";
    } else if(s_.squared_radius()==0){
      out << "(" << s_.center() << ")";
    } else {
      out << "(" << s_ << ", " << l_ << ": " 
	  << approximate_coordinate(Coordinate_index::X())
	  << " " << approximate_coordinate(Coordinate_index::Y())
	  << " " << approximate_coordinate(Coordinate_index::Z()) << ")";
    }
    return out;
  }



  CGAL::Comparison_result compare_c(const Sphere_line_intersection &o, Coordinate_index i) const;
  CGAL::Comparison_result compare_c(const Point_3 &o, Coordinate_index i) const;
  CGAL::Comparison_result compare_c( NT o, Coordinate_index i) const;

  const Quadratic_NT& exact_coordinate(Coordinate_index i) const {
    if (!has_exact_) {
      /*if (has_simple_coordinate(i)) {
	//std::cout << "Coordinate " << i << " is simple " << simple_coordinate(i) << std::endl;
	exact_[i.index()]= Quadratic_NT(simple_coordinate(i));
	} else {*/
 	set_has_exact(true);
	//} 
      //std::cout << "Coordinate " << i << " is " << exact_[i.index()] << std::endl;
    }
    return exact_[i.index()];
  }

  bool has_simple_coordinate(Coordinate_index i) const {
    return (l_.to_vector()[i.index()] ==0 || s_.squared_radius()==0);
  }


  const NT& simple_coordinate(Coordinate_index i) const {
    CGAL_assertion(has_simple_coordinate(i));
    if (s_.squared_radius()==0){
      CGAL_assertion(root_of_2(i)==s_.center()[i.index()]);
      return s_.center()[i.index()];
    } else {
      CGAL_assertion(root_of_2(i)==l_.point()[i.index()]);
      return l_.point()[i.index()];
    }
  }

  double approximate_coordinate(Coordinate_index i) const {
    // NOTE this sucks
    return CGAL::to_double(exact_coordinate(i));
  }

  double approximate_coordinate(int i) const {
    // NOTE this sucks
    return CGAL::to_double(exact_coordinate(Coordinate_index(i)));
  }

  std::pair<double, double> interval_coordinate(Coordinate_index i) const {
    return CGAL::to_interval(exact_coordinate(i));
  }

  std::pair<double, double> interval_coordinate(int i) const {
     return CGAL::to_interval(exact_coordinate(Coordinate_index(i)));
  }

  void swap(Sphere_line_intersection &o) {
    std::swap(s_, o.s_);
    std::swap(l_, o.l_);
  }

  Sphere_line_intersection flip_on(Coordinate_index i) const {
    Sphere_line_intersection r;
    CGAL_assertion(0);
    return r;
  }

  Sphere_3 s_;
  Line_3 l_;
  mutable bool has_exact_;
  mutable Quadratic_NT exact_[3];
};

CGAL_OUTPUT1(Sphere_line_intersection);




template <class K>
inline CGAL::Comparison_result Sphere_line_intersection<K>::compare_c(const Sphere_line_intersection<K> &o, Coordinate_index CC) const {
  CGAL_assertion(is_valid());
  CGAL_assertion(o.is_valid());
  CGAL_LOG_WRITE(Log::LOTS, write(LOG_STREAM) << std::endl);
  CGAL_LOG_WRITE(Log::LOTS, o.write(LOG_STREAM) << std::endl);
  //double dmc= CGAL::to_double(mc);
  //double doc= CGAL::to_double(oc);
  //std::cout << "Performing exact comparison " << mc << " vs " << oc << std::endl;
  Comparison_result ans= CGAL::compare(exact_coordinate(CC), o.exact_coordinate(CC));
#ifndef NDEBUG
  /*write(std::cout) << std::endl;
    o.write(std::cout) << std::endl;*/
 
  Comparison_result check= CGAL::compare(root_of_2(CC),
					 o.root_of_2(CC));
  if (ans != check) {
    CGAL_ERROR(ans << " " << check);
    CGAL_ERROR(exact_coordinate(CC) << ": " 
	       << Interval_nt<true>(to_interval(exact_coordinate(CC))));
    CGAL_ERROR(o.exact_coordinate(CC)
	       << Interval_nt<true>(to_interval(o.exact_coordinate(CC))));
    CGAL_ERROR(root_of_2(CC)
	       << Interval_nt<true>(to_interval(root_of_2(CC))));
    CGAL_ERROR(o.root_of_2(CC)
	       << Interval_nt<true>(to_interval(o.root_of_2(CC))));
    CGAL_assertion(0);
  }
#endif
  return ans;
  /*if (mc < oc) return CGAL::SMALLER;
  else if (oc < mc) return CGAL::LARGER;
  else return CGAL::EQUAL;*/
}

template <class K>
inline CGAL::Comparison_result Sphere_line_intersection<K>::compare_c(const Sphere_line_intersection<K>::Point_3 &o, Coordinate_index CC) const {
  CGAL_assertion(is_valid());
  //Quadratic_NT mc= exact_coordinate(CC);
  NT oc= o[CC.index()];
  return compare_c(oc, CC);
  //if (mc < oc) return CGAL::SMALLER;
  //else if (oc < mc) return CGAL::LARGER;
  //else return CGAL::EQUAL;
}


template <class K>
inline CGAL::Comparison_result Sphere_line_intersection<K>::compare_c(NT o, Coordinate_index CC) const {
  CGAL_assertion(is_valid());
  Quadratic_NT mc= exact_coordinate(CC);
  return CGAL::compare(mc, o);
  /*if (mc - o <0) return CGAL::SMALLER;
  else if (o < mc) return CGAL::LARGER;
  else return CGAL::EQUAL;*/
}
CGAL_AOS3_END_INTERNAL_NAMESPACE

#endif
