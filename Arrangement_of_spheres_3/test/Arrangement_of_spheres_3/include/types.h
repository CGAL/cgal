#ifndef TYPES_H_
#define TYPES_H_
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Gmpq NT;
typedef CGAL::Cartesian<NT> K;
typedef K::Sphere_3 Sphere;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Line_3 Line_3;
typedef K::Line_2 Line_2;
typedef K::Point_2 Point_2;
typedef K::Vector_2 Vector_2;
typedef K::Circle_2 Circle;
typedef K::Segment_2 Segment_2;

typedef CGAL::Cartesian<double> DK;
typedef DK::Point_2 DPoint;
typedef DK::Sphere_3 DSphere;

// lets be fancy bits are LRTBAI


struct Feature{
  /* for rules
     Inside means below a horizontal arc and to the left of a vertical one
  */
  enum PART_BITS {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, ARC_BIT=16, INF_BIT=32, IN_BIT=64};
  enum Part {INVALID=0, L_RULE=L_BIT, R_RULE=R_BIT, T_RULE=T_BIT, B_RULE=B_BIT,
	     LB_ARC=L_BIT|B_BIT|ARC_BIT, 
	     LT_ARC=L_BIT|T_BIT|ARC_BIT,
	     RT_ARC=R_BIT|T_BIT|ARC_BIT, 
	     RB_ARC=R_BIT|B_BIT|ARC_BIT, 
	     L_INF=L_BIT|INF_BIT,
	     R_INF=R_BIT|INF_BIT,
	     T_INF=T_BIT|INF_BIT, 
	     B_INF=B_BIT|INF_BIT};

  Feature(int i, Part pt): index_(i), pt_(pt){
    CGAL_precondition(is_finite());
  }
  Feature(Part pt): index_(-1), pt_(pt){
    CGAL_precondition(!is_finite());
  }
  Feature():index_(-1), pt_(INVALID){}
  bool is_valid() const {
    return pt_>0 && pt_ < 128
      //&& (pt_ ^ IN_BIT)
      && (index_>=0 || !is_finite()) 
      && !(is_arc() && !is_finite() );
  }
  Feature other_side() const {
    Feature ret=*this;
    ret.set_is_inside(!is_inside());
    CGAL_assertion(!(is_inside() && ret.is_inside()));
    return ret;
  }
  bool is_inside() const {
    return pt_&IN_BIT;
  }
  void set_is_inside(bool tf) {
    if (tf) pt_= pt_| IN_BIT;
    else pt_= pt_ & (~IN_BIT);
  }
  bool is_top() const {
    return pt_ & T_BIT;
  }
  bool is_left() const {
    return pt_ & L_BIT;
  }
  int index() const {
    return index_;
  }
  Part part() const {
    return static_cast<Part>(pt_);
  }
  bool is_rule() const {
    return ! is_arc();
  }
  bool is_negative() const {
    return pt_&L_BIT || pt_&B_BIT;
  }
  bool is_finite() const {
    return ! (pt_ & INF_BIT);
  }
  bool is_arc() const {
    return pt_ &ARC_BIT;
  }
  bool is_vertical() const {
    CGAL_precondition(is_rule());
    if (is_finite()) {
      return (pt_ & T_BIT) || (pt_ & B_BIT);
    } else {
      return (pt_ & L_BIT) || (pt_ & R_BIT);
    }
  }
  bool operator==(const Feature &o) const {
    return index_== o.index_ && pt_== o.pt_;
  }
  bool operator!=(const Feature &o) const {
    return index_!= o.index_ || pt_!= o.pt_;
  }
  bool operator<(const Feature &o) const {
    CGAL_precondition(is_valid() && o.is_valid());
    if (index_ < o.index_) return true;
    else if (index_ > o.index_) return false;
    else return pt_ < o.pt_;
  }
  bool operator>(const Feature &o) const {
    CGAL_precondition(is_valid() && o.is_valid());
    if (index_ > o.index_) return true;
    else if (index_ < o.index_) return false;
    else return pt_ > o.pt_;
  }
  std::ostream &write(std::ostream&out) const {
    CGAL_precondition(is_valid());
    out << to_string(pt_);
    if (is_finite())  out << " of " << index_;
    return out;
  }

  static const char *to_string(int pt){
    switch (pt) {
    case INVALID:
      return "Invalid";
    case L_RULE:
      return "L_rule";
    case R_RULE:
      return "R_rule";
    case T_RULE:
      return "T_rule";
    case B_RULE:
      return "B_rule";
    case LB_ARC:
      return "LB_arc";
    case LT_ARC:
      return "LT_arc";
    case RT_ARC:
      return "RT_arc"; 
    case RB_ARC:
      return "RB_arc"; 
    case L_INF:
      return "L_inf";
    case R_INF:
      return "R_inf";
    case T_INF:
      return "T_inf"; 
    case B_INF:
      return "B_inf";

    case L_RULE | IN_BIT :
      return "Li_rule";
    case R_RULE | IN_BIT:
      return "Ri_rule";
    case T_RULE | IN_BIT:
      return "Ti_rule";
    case B_RULE | IN_BIT:
      return "Bi_rule";
    case LB_ARC | IN_BIT:
      return "LBi_arc";
    case LT_ARC | IN_BIT:
      return "LTi_arc";
    case RT_ARC | IN_BIT:
      return "RTi_arc"; 
    case RB_ARC | IN_BIT:
      return "RBi_arc";
    case L_INF | IN_BIT:
      return "Li_inf";
    case R_INF | IN_BIT:
      return "Ri_inf";
    case T_INF | IN_BIT:
      return "Ti_inf"; 
    case B_INF | IN_BIT:
      return "Bi_inf";
    default:
      std::cerr << "Oops, I forgot: " << pt <<std::endl; 
      CGAL_assertion(0);
      return "Missing";
    }
  }
  const char *to_string() const {
    return to_string(pt_);
  }

private:
  int index_;
  int pt_;
};

inline std::ostream &operator<<(std::ostream &out, Feature f) {
  return f.write(out);
}

class Point{
public:
  enum Type {SS, SR, RR};
  Point(){}
  Point(Feature a, Feature b): a_(a), b_(b){
    //CGAL_precondition(a_.is_rule() || b_.is_rule());
    //a_= std::min(a,b);
    //b_=std::max(a,b);
  }
 

  bool operator==(const Point &o) const {
    return a_==o.a_ && b_==o.b_;
  }
  bool operator<(const Point &o) const {
    if (a_<o.a_) return true;
    else if (a_ > o.a_) return false;
    else return (b_ < o.b_);
  }

  std::ostream &write(std::ostream &out) const {
    out << "(" << a_ << ", " << b_ << ")";
    return out;
  }
  
  Type type() const {
    if (a_.is_rule() && b_.is_rule()) return RR;
    else if (a_.is_rule() || b_.is_rule()) return SR;
    else return SS;
  }

  Feature sphere(int i) const {
    if (i==0) {
      CGAL_precondition(type()== SR || type()== SS);
      if (!a_.is_rule()) return a_;
      else return b_;
    } else {
       CGAL_precondition( type()== SS);
       return b_;
    }
  }

  Feature rule(int i) const {
    if (i==0) {
      CGAL_precondition(type()== SR || type()== RR);
      if (a_.is_rule()) return a_;
      else return b_;
    } else {
       CGAL_precondition( type()== RR);
       return b_;
    }
  }

protected:
  Feature a_, b_;
};

inline std::ostream &operator<<(std::ostream &out, Point f) {
  return f.write(out);
}

struct Edge {
  
  Edge(Point s, Feature su, Point t): sup_(su), s_(s), t_(t){}

  bool operator<(const Edge &o) const {
    if (sup_ < o.sup_) return true;
    else if (o.sup_ < sup_) return false;
    else if (s_ < o.s_) return true;
    else if (o.s_ < s_) return false;
    else return t_ < o.t_;
  }

  Feature sup_;
  Point s_, t_;
};

#endif
