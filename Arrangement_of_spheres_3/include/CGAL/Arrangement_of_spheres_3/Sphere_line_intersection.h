#ifndef SPHERE_LINE_INTERSECTION_H
#define SPHERE_LINE_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/utilities.h>
#include <CGAL/tags.h>


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
  static typename K::Point_3 unproject(typename K::Point_2 p) {
    return typename K::Point_3(p.x(), p.y(),0);
  }
  static typename K::Vector_3 unproject(typename K::Vector_2 p) {
    return typename K::Vector_3(p.x(), p.y(),0);
  }
  static typename K::Line_3 unproject(typename K::Line_2 p) {
    return typename K::Line_3(unproject(p.point()), unproject(p.to_vector()));
  }
public:
  typedef K T;
  typedef typename K::FT NT;
  typedef typename CGAL::Root_of_traits<NT>::RootOf_2 Exact_NT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Line_3 Line_3;
  typedef typename K::Sphere_3 Sphere_3;

  //enum Type {FIRST, SECOND, ONE, INVALID};

  Sphere_line_intersection(): l_(Point_3(0,0,0), Vector_3(0,0,0)){
    CGAL_assertion(l_.is_degenerate());
  }

  /*Sphere_line_intersection(NT n):  type_(false){
    initialize_const(n);
    }*/
  Sphere_line_intersection(typename T::Point_3 p3): s_(p3,0), l_(p3, Vector_3(1,0,0)){
  }


  Sphere_line_intersection(typename T::Point_3 p3, typename T::Line_3 l): s_(p3,0),l_(l){
    CGAL_exactness_precondition(l.has_on(p3));
  }
  // assume the line hits
  Sphere_line_intersection( typename T::Sphere_3 s, 
			    typename T::Line_3 l, CGAL::Tag_true ): s_(s), 
								    l_(l){
  }
 
  Sphere_line_intersection(typename T::Sphere_3 s, typename T::Line_3 l): s_(s), l_(l) {
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

  /*Sphere_line_intersection(typename T::Line_3 l, bool pos): l_(l) {
    if (pos) {
      type_=PINF;
    } else {
      type_=NINF;
    }
    }*/

  // 2D constructors
  
  Sphere_line_intersection(typename T::Circle_2 s, typename T::Line_2 l): s_(unproject(s.center()),
					   s.squared_radius()), 
					l_(unproject(l)) {
    typename T::Point_2 cp= l.projection(s.center()); //typename T::closest_point(l, s.center());
    if (s.bounded_side(cp) == CGAL::ON_UNBOUNDED_SIDE){
      s_= Sphere_3();
      l_=Line_3(Point_3(0,0,0), Vector_3(0,0,0));
      CGAL_assertion(l_.is_degenerate());
    }
  }

  Sphere_line_intersection(typename T::Circle_2 s, typename T::Line_2 l, CGAL::Tag_true): s_(unproject(s.center()),
											     s.squared_radius()), 
											  l_(unproject(l)) {
  }

  /*Sphere_line_intersection(typename T::Line_2 l, bool pos): l_(unproject(l)) {
  if (pos) {
      type_=PINF;
    } else {
      type_=NINF;
    }
    }*/

  Sphere_line_intersection(typename T::Point_2 p3, typename T::Line_2 l): s_(unproject(p3),0),
									  l_(unproject(l)) {
  }

  /*bool is_finite() const {
    return type_ != NINF && type_!= PINF;
    }*/

  /*bool is_finite(int i) const {
    return is_finite() || l_.to_vector()[i]==0;
    }*/

  bool is_valid() const {
    return !l_.is_degenerate();
  }

  Sphere_3 sphere() const {
    return s_;
  }

  Line_3 line() const {
    return l_;
  }


  CGAL::Comparison_result compare_on_line(const Point_3 &pt) const {
    CGAL_precondition(line().has_on(pt));
    for (unsigned int i=0; i< 3; ++i){
      CGAL::Comparison_result c= compare(pt,i);
      if (c != CGAL::EQUAL) {
	if (line().to_vector()[i] >0) return c;
	else return CGAL::Comparison_result(-c);
      }
    }
    return CGAL::EQUAL;
  }

  /*void set_sphere(Sphere s) const {
    if (s_ != s){
      Sphere_line_intersection tmp1(s, l_, type_);
      if (!tmp1.is_valid()) return false;
      else {
	return *this == tmp1;
      }
    } else {
      return true;
    }
    }*/

  std::ostream &write(std::ostream &out) const {
    //out << "(" << s_ << ", " << l_ << ", (" << lb_ << "..." << ub_ << "))";
    if (!is_valid()) {
      out << "INVALID";
    } else if(s_.squared_radius()==0){
      out << "(" << s_.center() << ")";
    } else {
      out << "(" << s_ << ", " << l_ << ": " 
	  << approximate_coordinate(0)
	  << " " << approximate_coordinate(1)
	  << " " << approximate_coordinate(2) << ")";
    }
    return out;
  }

  /*bool is_equivalent_sphere(Sphere s) const {
    std::cout<< "orig " << s_ << " asking about " << s << std::endl;
    Sphere_line_intersection tmp0=*this, tmp1=*this;
    tmp0.s_=s;
    if ( tmp0 != tmp1) {
      std::cout << tmp0.double_value() << " vs " << tmp1.double_value() << std::endl;
      return false;
    } else return true;
    }*/

  CGAL::Comparison_result compare(const Sphere_line_intersection &o, int CC) const;
  CGAL::Comparison_result compare(const Point_3 &o, int CC) const;

  Exact_NT exact_coordinate(int i) const {
    //CGAL_precondition(is_finite(i));
    if (l_.to_vector()[i]==0) return Exact_NT(l_.point()[i]);
    else if (s_.squared_radius()==0) return s_.center()[i];
    else return exact_intersection<K>(s_, l_, i);
  }

 bool has_simple_coordinate(int i){
    return (l_.to_vector()[i] ==0 || s_.squared_radius()==0);
  }


  NT simple_coordinate(int i){
    CGAL_assertion(has_simple_coordinate(i));
    if (s_.squared_radius()==0) return s_.center()[i];
    else return l_.point()[i];
  }

  /*Exact_NT clipped_exact_coordinate(int i, NT max=NT(1000)) const {
    CGAL_precondition(max>0);
    if (is_finite(i)) return exact_coordinate(i);
    else if (l_.to_vector()[i]==0) return Exact_NT(l_.point()[i]);
    else if (type_==NINF) return Exact_NT(-max);
    else return Exact_NT(max);
    }*/

  double approximate_coordinate(int i) const {
    /*if (!is_finite(i)) {
      if (type_== PINF) return std::numeric_limits<double>::infinity();
      else return -std::numeric_limits<double>::infinity();
      } else {*/
    return CGAL::to_double(exact_coordinate(i));
      //}
  }

  void swap(Sphere_line_intersection &o) {
    std::swap(s_, o.s_);
    std::swap(l_, o.l_);
  }

  Sphere_line_intersection flip_on(int c) const {
    Sphere_line_intersection r;
    CGAL_assertion(0);
    return r;
  }

  Sphere_3 s_;
  Line_3 l_;
};

template <class K>
inline std::ostream &operator<<(std::ostream &out, const Sphere_line_intersection<K> &r){
  return r.write(out);
}



template <class K>
inline CGAL::Comparison_result Sphere_line_intersection<K>::compare(const Sphere_line_intersection<K> &o, int CC) const {
  CGAL_assertion(is_valid());
  CGAL_assertion(o.is_valid());
  Exact_NT mc= exact_coordinate(CC);
  Exact_NT oc= o.exact_coordinate(CC);
  //std::cout << "Performing exact comparison " << mc << " vs " << oc << std::endl;
  if (mc < oc) return CGAL::SMALLER;
  else if (oc < mc) return CGAL::LARGER;
  else return CGAL::EQUAL;
}

template <class K>
inline CGAL::Comparison_result Sphere_line_intersection<K>::compare(const Sphere_line_intersection<K>::Point_3 &o, int CC) const {
  CGAL_assertion(is_valid());
  Exact_NT mc= exact_coordinate(CC);
  NT oc= o[CC];
  if (mc < oc) return CGAL::SMALLER;
  else if (oc < mc) return CGAL::LARGER;
  else return CGAL::EQUAL;
}


#endif
