#ifndef ARRANGEMENT_SPHERE_3_VERTEX_H
#define ARRANGEMENT_SPHERE_3_VERTEX_H

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>


/*
  If the two curves are both arcs, then the vertex is the one that is
  ccw with respect to the ordered centers.

  If they are both rules then the vertical one is first.
  If there is an arc and a rule then the rule is first.
*/
class Combinatorial_vertex{
public:
  enum Type {SS, SR, RR, SPECIAL};
  
  static Combinatorial_vertex make_special(Combinatorial_curve::Key i) {
    Combinatorial_vertex cv;
    cv.a_= Combinatorial_curve::make_special(i);
    cv.b_= Combinatorial_curve::make_special(i);
    return cv;
  }

  Combinatorial_vertex(){}
  Combinatorial_vertex(Combinatorial_curve a, Combinatorial_curve b){
    if (a.is_rule() && b.is_rule() && a.is_vertical()
	|| b.is_arc()) {
      a_=a.strip_inside();
      b_=b.strip_inside();
    } else {
      b_=a.strip_inside();
      a_=b.strip_inside();
    }

    CGAL_postcondition(is_valid());
    //CGAL_precondition(a_.is_rule() || b_.is_rule());
    //a_= std::min(a,b);
    //b_=std::max(a,b);
  }
 

  bool is_valid() const {
    if (a_.is_rule() && b_.is_rule() && !a_.is_vertical()) return false;
    if (b_.is_rule() && !a_.is_rule()) return false;
    return true;
  }

  bool operator==(const Combinatorial_vertex &o) const {
    return a_==o.a_ && b_==o.b_;
  }
  bool operator!=(const Combinatorial_vertex &o) const {
    return a_!=o.a_ || b_!=o.b_;
  }
  bool operator<(const Combinatorial_vertex &o) const {
    if (a_<o.a_) return true;
    else if (a_ > o.a_) return false;
    else return (b_ < o.b_);
  }

  std::ostream &write(std::ostream &out) const {
    out << "(" << a_ << ", " << b_ << ")";
    return out;
  }
  bool is_special() const {
    return a_.is_special() && b_.is_special();
  }
  Type type() const {
    if (b_.is_rule()) return RR;
    else if (a_.is_rule()) return SR;
    else return SS;
  }

  Combinatorial_curve sphere(int i) const {
    CGAL_precondition(type()== SR || type()== SS);
    if (i==0 && type() == SS) {
      return a_;
    } else {
      return b_;
    }
  }

  void replace_rule(Combinatorial_curve c) {
    if (type() == SR) {
      a_= c;
    } else if (type()== RR) {
      if (c.is_vertical() && a_.is_vertical()){
	a_=c;
      } else {
	b_=c;
      }
    } else {
      CGAL_assertion(0);
    }
  }

  Combinatorial_curve rule(int i) const {
    if (i==0) {
      return a_;
    } else {
      CGAL_precondition( type()== RR);
      return b_;
    }
  }

  Combinatorial_curve first() const {
    return a_;
  }

  Combinatorial_curve second() const {
    return b_;
  }
  bool is_finite() const {
    return a_.is_finite() && b_.is_finite();
  }

protected:
  Combinatorial_curve a_, b_;
};

inline std::ostream &operator<<(std::ostream &out, Combinatorial_vertex f) {
  return f.write(out);
}

#endif
