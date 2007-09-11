#ifndef FILTERED_SPHERE_LINE_INTERSECTION_H
#define FILTERED_SPHERE_LINE_INTERSECTION_H


#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <iostream>
//#include <CGAL/Arrangement_of_spheres_3/utilities.h>



//using doubles for ease of implementation, may not be the fastest
/*
  cases to handle:
  - tangent line -- fall back to algebraic since interval is small
  - tiny sphere--fall back to algebraic since interval is small
  - missing line
  - vertical line-- asserted
   
  I want to make sure that the interval is open. Now it is not. 
*/



CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


//int num_exact_tests_=0;
//int num_tests_=0;


template <class K, class C>
struct Filtered_sphere_line_intersection: public Sphere_line_intersection<K> {
  enum Extended_comparison_result {UNKNOWN = -2, SMALLER= CGAL::SMALLER, 
				   EQUAL= CGAL::EQUAL, LARGER= CGAL::LARGER};
  typedef Sphere_line_intersection<K> P;
  typedef Filtered_sphere_line_intersection<K,C> This;
  


  static typename P::Coordinate_index coordinate() {
    //return Coordinate_index(C);
    switch (C::index()) {
    case 0: return Coordinate_index::X();
    case 1: return Coordinate_index::Y();
    default: return Coordinate_index::Z();
    }
  }


  Filtered_sphere_line_intersection(): lb_(1), ub_(-1){
  }


  Filtered_sphere_line_intersection(typename K::Point_3 p3,
				    typename K::Line_3 l): 
    P(p3, l){
    std::pair<double,double> i= CGAL::to_interval(p3[C::index()]);
    lb_=i.first;
    ub_=i.second;
  }


 
 
  Filtered_sphere_line_intersection(typename K::Sphere_3 s, 
				    typename K::Line_3 l): P(s,l){
    // hope the compiler collapses the duplicate work
    if (!P::is_valid()){
      lb_=std::numeric_limits<double>::infinity();
      ub_=lb_;
    } else {
      typename K::Point_3 cp= l.projection(s.center()); //K::closest_point(l, s.center());
      initialize(cp[C::index()]);
    }
  }

 
  Filtered_sphere_line_intersection(typename P::NT p): P(typename K::Point_3(p,p,p),
							 typename K::Line_3(typename K::Point_3(p,p,p), 
									    typename K::Vector_3(1,1,1))){
    std::pair<double,double> i= CGAL::to_interval(p);
    lb_=i.first;
    ub_=i.second;
  }

  //Sphere_line_intersection(double lb, double ub):lb_(lb), ub_(ub){}

  CGAL_GETNR(double, lb, return lb_);
  CGAL_GETNR(double, ub, return ub_);


  static Filtered_sphere_line_intersection infinity_rep() {
    typename K::Point_3 p(0,0,0);
    typename P::NT vs[]={0,0,0};
    vs[C::index()]=1;
    typename P::Vector_3 v(vs[0], vs[1], vs[2]);
    typename P::Line_3 l(p,v);
    return Filtered_sphere_line_intersection(l, true);
  }


  CGAL_COMPARISONS_COMPARE;

  void initialize() {
    // fix this
    typename P::NT t= (P::line().to_vector()*(P::sphere().center()-CGAL::ORIGIN) 
		       - P::line().to_vector()*(P::line().point()-CGAL::ORIGIN))/(P::line().to_vector()*P::line().to_vector());
    typename P::NT r= P::line().point()[C::index()]+t*P::line().to_vector()[C::index()];
    initialize(r);
  }

  void initialize(typename P::NT cc) {
    //NT cc =closest_coord<C>(s_.center(), P::line());
    /*CGAL_precondition_code(typename P::NT t=(cc-P::line().point()[C])/P::line().to_vector()[C]);
    CGAL_precondition_code(typename P::Point_3 pol=P::line().point()+t*P::line().to_vector());
    CGAL_precondition(P::sphere().bounded_side(pol) != CGAL::ON_UNBOUNDED_SIDE);*/
    
    
    std::pair<double,double> i=CGAL::to_interval(cc);
    if (P::line().to_vector()[C::index()]>0) {
      ub_=i.second;
      lb_=CGAL::to_interval(P::sphere().center()[C::index()] - (std::max)(P::sphere().squared_radius(), typename P::NT(1))).first;
    } else {
      lb_=i.first;
      ub_=CGAL::to_interval(P::sphere().center()[C::index()] + (std::max)(P::sphere().squared_radius(), typename P::NT(1))).second;
    }
  }

  /*void initialize_const(NT n) {
    std::pair<double, double> i= CGAL::to_interval(n);
    lb_=i.first;
    ub_=i.second;
    if (lb_ != ub_) {
    // initialize sphere and line
    NT v[3]={0,0,0};
    v[C::index()]=1;
    l_=Line_3(Point_3(0,0,0),Vector_3(v[0], v[1], v[2]));
    NT r= n-NT(ub_);
    v[C::index()]=n+r;
    s_=Sphere(Point_3(v[0], v[1], v[2]), r*r);
    }
    }*/

  std::ostream &write(std::ostream &out) const {
    //out << "(" << s_ << ", " << l_ << ", (" << lb_ << "..." << ub_ << "))";
    P::write(out) << " [" << lb_ << "..." << ub_ << "] = " << double_value();
    return out;
  }

  double double_value() const {
    if (ub_== lb_) return ub_;
    else if (ub_-lb_ < .00001) return .5*(lb_+ub_);
    else {
      double d=CGAL::to_double(P::exact_coordinate(coordinate()));
      CGAL_assertion(d>=lb_ && d<= ub_);
      return d;
    }
  }

  std::pair<double, double> interval_value() const {
    return std::make_pair(lb_, ub_);
  }

  void refine(double bm) const {
    if (bm==lb_ || bm==ub_) return;
    CGAL_precondition( bm > lb_ && bm < ub_);
    if (P::line().to_vector()[C::index()] == 0){
      if (typename P::NT(bm) < P::line().point()[C::index()]){
	lb_=bm;
      } else {
	ub_=bm;
      }
    } else {
      typename P::NT t=(typename P::NT(bm)-P::line().point()[C::index()])
	/P::line().to_vector()[C::index()];
      typename K::Point_3 p=  P::line().point()+t*P::line().to_vector();

      //= point_on_line<K, C>(P::line(), typename P::NT(bm));
      CGAL::Bounded_side bs= P::sphere().bounded_side(p);
      switch(bs) {
      case CGAL::ON_BOUNDED_SIDE: if (P::line().to_vector()[C::index()] > 0) ub_=bm; else lb_=bm; break;
      case CGAL::ON_UNBOUNDED_SIDE: if (P::line().to_vector()[C::index()] > 0) lb_=bm; else ub_=bm; break;
      case CGAL::ON_BOUNDARY: lb_=ub_=bm; break;
      };
    }
    /*{
      CGAL_assertion(P::exact_coordinate(typename P::Coordinate_index(C)) 
		     >= typename P::Quadratic_NT(typename P::NT(lb_)));
      CGAL_assertion(P::exact_coordinate(typename P::Coordinate_index(C)) 
		     <= typename P::Quadratic_NT(typename P::NT(ub_)));
		     }*/
  }

  

  void refine() const {
    refine( .5*(lb_+ub_));
  }

  This operator-() const {
    CGAL_assertion(0);
    return This();
  }

  static void test();

  
  typename P::NT coord_on_line(typename P::NT v,
			       typename P::Coordinate_index CC) const {
    typename P::NT t=(v-P::line().point()[C::index()])/P::line().to_vector()[C::index()];
    return  P::line().point()[CC.index()]+t*P::line().to_vector()[CC.index()];
  }

  /*void assert_is_valid() {
    if (lb_ != ub_ && .5*(lb_+ub_) != lb_ && .5*(lb_+ub_) != ub_) {
    if (first_) {
    CGAL_assertion(bounded_side_of_sphere(s_, point_on_line<C>(NT(ub_)))!= CGAL::ON_UNBOUNDED_SIDE);
    CGAL_assertion(bounded_side_of_sphere(s_, point_on_line<C>(NT(lb_)))!= CGAL::ON_BOUNDED_SIDE);
    } else {
    CGAL_assertion(bounded_side_of_sphere(s_, point_on_line<C>(NT(ub_)))!= CGAL::ON_BOUNDED_SIDE);
    CGAL_assertion(bounded_side_of_sphere(s_, point_on_line<C>(NT(lb_)))!= CGAL::ON_UNBOUNDED_SIDE);
    }
    }
    }*/

  void swap(This &o) {
    std::swap(lb_, o.lb_);
    std::swap(ub_, o.ub_);
    P::swap(o);
  }


  static Extended_comparison_result compare_bounds(const This &a,
						   const This &b,
						   typename P::Coordinate_index CC){
    //typedef Filtered_sphere_line_intersection<K,C> FS;
    if (a.line().to_vector()[CC.index()] == 0 || b.line().to_vector()[CC.index()]== 0){
      // just punt and fall back on exact for now. 
      return UNKNOWN;
    } else {
      typename P::NT tl= a.coord_on_line(typename P::NT(a.lb()), CC);
      typename P::NT tu= a.coord_on_line(typename P::NT(a.ub()), CC);
      typename P::NT ol= b.coord_on_line(typename P::NT(b.lb()), CC);
      typename P::NT ou= b.coord_on_line(typename P::NT(b.ub()), CC);
      if (tl > tu) std::swap(tl,tu);
      if (ol > ou) std::swap(ol,ou);
      if (tu < ol) {
	CGAL_assertion(C::object() != CC 
		       || a.double_value() < b.double_value());
	return SMALLER;
      } else if (ou < tl) {
	CGAL_assertion(C::object() != CC 
		       || a.double_value() > b.double_value());
	return LARGER;
      }
      else if (tl==tu && ol==ou && tl==ol) return EQUAL;
      else return UNKNOWN; 
    }
  }


  static Extended_comparison_result compare_bounds_C(const This  &a,
						     const This &b){
    if (a.ub() < b.lb()) {
      CGAL_assertion(a.double_value() < b.double_value());
      return SMALLER;
    }
    else if (a.lb() > b.ub()){
      CGAL_assertion(a.double_value() > b.double_value());
      return LARGER;
    } else if (a.lb()== b.lb() && a.lb()==a.ub() && b.lb() == b.ub()) return EQUAL;
    else return UNKNOWN;
  }
  

  inline CGAL::Comparison_result compare_2(const This &o,
					   typename P::Coordinate_index CC, 
					   int sct=0) const {
    Extended_comparison_result ecr= compare_bounds(*this, o, CC);
    if (ecr != UNKNOWN) return CGAL::enum_cast<CGAL::Comparison_result>(ecr); 
    else {
      if (sct < 5) {
	refine();
	o.refine();
	return compare_2(o, CC, sct+1);
      } else {
	//++num_exact_tests_;
	return P::compare_c(o, CC);
      }
    }
  }


  inline CGAL::Comparison_result compare_C2(const This &o, int sct=0) const {
    CGAL_precondition(P::is_valid() && o.is_valid());
    Extended_comparison_result ecr= compare_bounds_C(*this, o);
    if (ecr != UNKNOWN) return CGAL::enum_cast<CGAL::Comparison_result>(ecr); 
    else {
      if (sct < 5) {
	refine();
	o.refine();
	return compare_C2(o, sct+1);
      } else {
	//++num_exact_tests_;
	return static_cast<const P*>(this)->compare(o, coordinate());
      }
    }
  }
  inline CGAL::Comparison_result  compare(const This &o, 
					  typename P::Coordinate_index CC=typename P::Coordinate_index(Coordinate_Z())) const {
    //++num_tests_;
    // see if intervals overlap
    Extended_comparison_result ecr= compare_bounds(*this, o, CC);
    if (ecr != UNKNOWN) return CGAL::enum_cast<CGAL::Comparison_result>(ecr);
    // question of <= vs <
    // if so cut at shared part
    // not useful if we are comparing at another coordinate than C, but not really harmful either.
    else if (C::object()==CC && lb_ <= o.lb_) {
      if (ub_ != o.lb_) refine(o.lb_);
      if (ub_ < o.ub_) o.refine(ub_);
    } else if (C::object()==CC && o.lb_ <= lb_) {
      if (o.ub_ != lb_) o.refine(lb_);
      if (o.ub_ < ub_) refine(o.ub_);
    }
    return compare_2(o, CC);
    // if still, subdivide a few times
    // finally use algebraic numbers
  }

  
  inline CGAL::Comparison_result  compare_C(const This &o) const {
    //++num_tests_;
    // see if intervals overlap
    Extended_comparison_result ecr= compare_bounds_C(*this, o);
    if (ecr != UNKNOWN) return CGAL::enum_cast<CGAL::Comparison_result>(ecr);
    else if (lb_ <= o.lb_) {
      if (ub_ != o.lb_) refine(o.lb_);
      if (ub_ < o.ub_) o.refine(ub_);
    } else if (o.lb_ <= lb_) {
      if (o.ub_ != lb_) o.refine(lb_);
      if (o.ub_ < ub_) refine(o.ub_);
    }
    return compare_C2(o);
  }


  mutable double lb_, ub_;
};


CGAL_OUTPUT2(Filtered_sphere_line_intersection);
CGAL_SWAP2(Filtered_sphere_line_intersection);

CGAL_AOS3_END_INTERNAL_NAMESPACE

CGAL_REAL_EMBEDDABLE2(CGAL_AOS3_INTERNAL_NS::Filtered_sphere_line_intersection);
CGAL_HAS_INFINITY2(CGAL_AOS3_INTERNAL_NS::Filtered_sphere_line_intersection);





#endif
