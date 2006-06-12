#ifndef CGAL_ARRANGEMENT_SPHERES_3_FEATURE_H
#define CGAL_ARRANGEMENT_SPHERES_3_FEATURE_H

#include <CGAL/Arrangement_of_spheres_traits_3.h>

struct Combinatorial_curve{
  /* for rules
     Inside means below a horizontal arc and to the left of a vertical one
  */
  enum PART_BITS {L_BIT=1, R_BIT=2, T_BIT=4, B_BIT=8, ARC_BIT=16, /*INF_BIT=32,*/ IN_BIT=32};
  enum LOCATION_BITS {lL_BIT=1, lR_BIT=2, lT_BIT=4, lB_BIT=8, lIN_BIT=16, lOUT_BIT=32};
  enum Part {INVALID=0, L_RULE=L_BIT, R_RULE=R_BIT, T_RULE=T_BIT, B_RULE=B_BIT,
	     LB_ARC=L_BIT|B_BIT|ARC_BIT, 
	     LT_ARC=L_BIT|T_BIT|ARC_BIT,
	     RT_ARC=R_BIT|T_BIT|ARC_BIT, 
	     RB_ARC=R_BIT|B_BIT|ARC_BIT/*, 
	     L_INF=L_BIT|INF_BIT,
	     R_INF=R_BIT|INF_BIT,
	     T_INF=T_BIT|INF_BIT, 
	     B_INF=B_BIT|INF_BIT*/};

  typedef Arrangement_of_spheres_traits_3::Key Key;

  Combinatorial_curve(int i, Part pt): index_(i), pt_(pt){
    CGAL_precondition(i>=0);
    CGAL_precondition(is_finite());
  }
  Combinatorial_curve(Key i, Part pt): index_(i), pt_(pt){
    //CGAL_precondition(is_finite());
  }
  /*explicit Combinatorial_curve(Part pt): pt_(pt){
    if (pt_== T_BIT || pt_ == R_BIT) {
      index_= Key(Key::TR);
    } else {
      CGAL_assertion(pt_ == B_BIT || pt_ == L_BIT);
      index_= Key(Key::BL);
    }
    CGAL_precondition(!is_finite());
    }*/
  Combinatorial_curve():pt_(INVALID){}
  bool is_valid() const {
    return pt_>0 && pt_ < 128
      //&& (pt_ ^ IN_BIT)
      && index_.is_valid() 
      && !(is_arc() && !is_finite() );
  }
  Combinatorial_curve other_side() const {
    Combinatorial_curve ret=*this;
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

  Combinatorial_curve strip_inside() const {
    Combinatorial_curve ret(index_, pt_& (~IN_BIT), true);
    CGAL_assertion(!ret.is_inside());
    CGAL_assertion(ret.is_valid());
    return ret;
  }
  bool is_top() const {
    return pt_ & T_BIT;
  }
  bool is_right() const {
    return pt_ & R_BIT;
  }
  bool is_bottom() const {
    return pt_ & B_BIT;
  }
  bool is_left() const {
    return pt_ & L_BIT;
  }
  /*int index() const {
    return index_.to_index();
    }*/
  Key key() const {
    return index_;
  }
  int constant_coordinate() const {
    CGAL_precondition(is_rule());
    int r;
    if ((pt_&R_BIT) || (pt_ &L_BIT)) r= 1;
    else r= 0;
    //if (!is_finite()) r=1-r;
    return r;
  }
  bool can_intersect(int t) const {
    CGAL_precondition( (t & ~(L_BIT | R_BIT | T_BIT | B_BIT)) == 0);
    if (is_rule()) {
      if (!is_vertical() && ((t&R_BIT) || (t &L_BIT))) return false;
      if (is_vertical() && ((t&T_BIT) || (t &B_BIT))) return false;
      if (((t & R_BIT) || (t & T_BIT)) && !is_inside()) return false;
      if (((t & L_BIT) || (t & B_BIT)) && is_inside()) return false;
      return true;
    } else {
      // can simplify
      if (t&R_BIT) return (is_right() && is_inside() || !is_right() && !is_inside());
      if (t&L_BIT) return (is_right() && !is_inside() || !is_right() && is_inside());
      if (t&T_BIT) return (is_top() && is_inside() || !is_top() && !is_inside());
      if (t&B_BIT) return (is_top() && is_inside() || !is_top() && !is_inside());
      CGAL_assertion(0);
      return false;
    }
  }

  int quadrant() const {
    if (is_rule()) {
      if (is_vertical() && is_inside()) return R_BIT;
      else if (is_vertical() && !is_inside()) return L_BIT;
      else if (!is_vertical() && is_inside()) return T_BIT;
      else return B_BIT;
    } else {
      if (is_inside()) {
	return pt_ & (T_BIT | B_BIT | L_BIT | R_BIT);
      } else {
	return (~pt_) & (T_BIT | B_BIT | L_BIT | R_BIT);
      }
    }
  }

  Part part() const {
    bool dont_use_part;
    return static_cast<Part>(pt_);
  }


  bool is_rule() const {
    return ! is_arc();
  }
  bool is_negative() const {
    return pt_&L_BIT || pt_&B_BIT;
  }
  bool is_finite() const {
    return index_.is_input(); //! (pt_ & INF_BIT);
  }
  bool is_arc() const {
    bool nb= (pt_ &ARC_BIT) != 0;
    if (pt_&ARC_BIT) {
      CGAL_assertion(nb);
    }
    if (!(pt_&ARC_BIT)) {
      CGAL_assertion(!nb);
    }
    return nb;
  }
  bool is_vertical() const {
    CGAL_precondition(is_rule());
    //if (is_finite()) {
      return (pt_ & T_BIT) || (pt_ & B_BIT);
      /*} else {
      return (pt_ & L_BIT) || (pt_ & R_BIT);
      }*/
  }
  bool is_same_side(Combinatorial_curve o) const {
    int u= o.pt_ & pt_; 
    return u &( T_BIT | L_BIT | R_BIT | B_BIT);
  }
  
  bool operator==(const Combinatorial_curve &o) const {
    return index_== o.index_ && pt_== o.pt_;
  }
  bool operator!=(const Combinatorial_curve &o) const {
    return index_!= o.index_ || pt_!= o.pt_;
  }
  bool operator<(const Combinatorial_curve &o) const {
    CGAL_precondition(is_valid() && o.is_valid());
    if (index_ < o.index_) return true;
    else if (index_ > o.index_) return false;
    else return pt_ < o.pt_;
  }
  bool operator>(const Combinatorial_curve &o) const {
    CGAL_precondition(is_valid() && o.is_valid());
    if (index_ > o.index_) return true;
    else if (index_ < o.index_) return false;
    else return pt_ > o.pt_;
  }
  std::ostream &write(std::ostream&out) const {
    CGAL_precondition(is_valid());
    out << to_string(pt_);
    /*if (is_finite())*/  out << index_;
    return out;
  }

  bool is_compatible_location(int i) const {
    //typedef Arrangement_of_spheres_traits_3::Sphere_location SL;
    if (is_rule()) {
      if (is_vertical()) {
	if ( !(i & lL_BIT) && is_inside()) return false;
	else if (  !(i & lR_BIT) && !is_inside()) return false;
	else return true;
      } else {
	if (  !(i & lB_BIT) && is_inside()) return false;
	else if ( !(i & lT_BIT) && !is_inside()) return false;
	else return true;
      }
    } else {
      if (is_inside()) {
	if (!(i & lIN_BIT)) return false;
	else return true;
      } else {
	if (!(i & lOUT_BIT)) return false;
	if (!(i & pt_)) return false;
	//CGAL_assertion(static_cast<int>(R_BIT) == static_cast<int>(R_BIT));
	//CGAL_assertion(static_cast<int>(L_BIT) == static_cast<int>(L_BIT));
	//CGAL_assertion(static_cast<int>(T_BIT) == static_cast<int>(T_BIT));
	//CGAL_assertion(static_cast<int>(B_BIT) == static_cast<int>(B_BIT));
	return true;
      }
    }
  }

  bool can_intersect(const Combinatorial_curve &o) const {
    CGAL_assertion(is_rule());
    return pt_&o.quadrant();
  }

  int is_weakly_incompatible(int i) const {
    int a= i&pt_;
    if (a== L_BIT || a==R_BIT) return 1;
    else if (a== T_BIT || a== B_BIT) return 0;
    else return -1;
  }

  static const char *to_string(int pt){
    switch (pt) {
    case INVALID:
      return "Invalid";
    case L_RULE:
      return "L";
    case R_RULE:
      return "R";
    case T_RULE:
      return "T";
    case B_RULE:
      return "B";
    case LB_ARC:
      return "LB";
    case LT_ARC:
      return "LT";
    case RT_ARC:
      return "RT"; 
    case RB_ARC:
      return "RB"; 
      /*case L_INF:
      return "L_inf";
    case R_INF:
      return "R_inf";
    case T_INF:
      return "T_inf"; 
    case B_INF:
    return "B_inf";*/

    case L_RULE | IN_BIT :
      return "Li";
    case R_RULE | IN_BIT:
      return "Ri";
    case T_RULE | IN_BIT:
      return "Ti";
    case B_RULE | IN_BIT:
      return "Bi";
    case LB_ARC | IN_BIT:
      return "LBi";
    case LT_ARC | IN_BIT:
      return "LTi";
    case RT_ARC | IN_BIT:
      return "RTi"; 
    case RB_ARC | IN_BIT:
      return "RBi";
      /*case L_INF | IN_BIT:
      return "Li_inf";
    case R_INF | IN_BIT:
      return "Ri_inf";
    case T_INF | IN_BIT:
      return "Ti_inf"; 
    case B_INF | IN_BIT:
    return "Bi_inf";*/
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

  Combinatorial_curve(Key i, int pt, bool): index_(i), pt_(pt){
  }

  Key index_;
  int pt_;
};

inline std::ostream &operator<<(std::ostream &out, Combinatorial_curve f) {
  return f.write(out);
}

#endif
