#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

#define inline 

inline bool Combinatorial_curve::can_intersect(int t) const {
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


inline int Combinatorial_curve::quadrant() const {
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


inline Rule_direction Combinatorial_curve::rule_direction(const Combinatorial_curve &a,
							  const Combinatorial_curve &b) {
  CGAL_precondition(a.is_arc());
  CGAL_precondition(b.is_arc());
  CGAL_precondition(a.key() == b.key());
  /*switch(a.pt_ & b.pt_ & (L_BIT | T_BIT | R_BIT | B_BIT)) {
    case R_BIT:
    return 0;
    case T_BIT:
    return 1;
    case L_BIT:
    return 2;
    case B_BIT:
    return 3;
    default:
    CGAL_assertion(0);
    return -1;
    }*/
  return Rule_direction::make_from_part(a.pt_ & b.pt_ 
					& (L_BIT | T_BIT | R_BIT | B_BIT));
}

  

inline bool Combinatorial_curve::is_compatible_location(int i) const {
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

inline const char *Combinatorial_curve::to_string(int pt){
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
  case SPECIAL:
    return "SPECIAL";
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


inline bool Combinatorial_curve::is_arc() const {
  bool nb= (pt_ &ARC_BIT) != 0;
  if (pt_&ARC_BIT) {
    CGAL_assertion(nb);
  }
  if (!(pt_&ARC_BIT)) {
    CGAL_assertion(!nb);
  }
  return nb;
}
inline bool Combinatorial_curve::is_vertical() const {
  CGAL_precondition(is_rule());
  //if (is_finite()) {
  return (pt_ & T_BIT) || (pt_ & B_BIT);
  /*} else {
    return (pt_ & L_BIT) || (pt_ & R_BIT);
    }*/
}
inline bool Combinatorial_curve::is_same_side(Combinatorial_curve o) const {
  int u= o.pt_ & pt_; 
  return u &( T_BIT | L_BIT | R_BIT | B_BIT);
}


inline std::ostream &Combinatorial_curve::write(std::ostream&out) const {
  //CGAL_precondition(is_valid());
  out << to_string(pt_);
  /*if (is_finite())*/  out << index_;
  return out;
}



inline Coordinate_index Combinatorial_curve::is_weakly_incompatible(int i) const {
  int a= i&pt_&(~lOUT_BIT) &(~lIN_BIT);
  if (a == L_BIT || a == R_BIT) return plane_coordinate(1);
  else if (a == T_BIT || a == B_BIT) return plane_coordinate(0);
  else return Coordinate_index();
}

inline Rule_direction Combinatorial_curve::rule_direction() const {
  CGAL_precondition(is_rule());
  return Rule_direction::make_from_part(pt_&(R_BIT | L_BIT | T_BIT|B_BIT));
}

inline int Combinatorial_curve::arc_index() const {
  CGAL_precondition(is_arc());
  if ((pt_& RT_ARC) == RT_ARC) return 0;
  else if ((pt_& LT_ARC) == LT_ARC) return 1;
  else if ((pt_& LB_ARC) == LB_ARC) return 2;
  else if ((pt_& RB_ARC) == RB_ARC) return 3;
  else {
    CGAL_assertion(0);
    return -1;
  }
}

inline Combinatorial_curve Combinatorial_curve::make_rule(Key k,
							  Rule_direction ruleindex) {
  return Combinatorial_curve(k, Part(ruleindex.part()));
}

inline Coordinate_index
Combinatorial_curve::constant_coordinate() const {
  CGAL_precondition(is_rule());
  Coordinate_index r;
  if ((pt_&R_BIT) || (pt_ &L_BIT)) 
    r= plane_coordinate(1);
  else r= plane_coordinate(0);
  //if (!is_finite()) r=1-r;
  //std::cout << "Constant coordinate of " << *this << " is " << r << std::endl;
  return r;
}

inline bool Combinatorial_curve::is_same_part(const Combinatorial_curve &o) const {
  CGAL_precondition(is_inside() == o.is_inside());
  CGAL_precondition(is_rule());
  CGAL_precondition(o.is_rule());
  return o.pt_==pt_;
}

inline void Combinatorial_curve::flip_rule(Key k) {
  CGAL_precondition(is_rule());
  index_=k;
  if (is_vertical()) {
    pt_ = pt_^ (T_BIT | B_BIT);
  } else {
    pt_ = pt_^ (L_BIT | R_BIT);
  }
}


inline const Combinatorial_curve Combinatorial_curve::strip_inside() const {
  Combinatorial_curve ret(index_, pt_& (~IN_BIT), true);
  CGAL_assertion(!ret.is_inside());
  CGAL_assertion(ret.is_valid());
  return ret;
}

inline const Combinatorial_curve Combinatorial_curve::other_side() const {
  Combinatorial_curve ret=*this;
  ret.set_is_inside(!is_inside());
  CGAL_assertion(!(is_inside() && ret.is_inside()));
  return ret;
}


inline bool Combinatorial_curve::is_valid() const {
  return pt_>0 && pt_ < 128
    //&& (pt_ ^ IN_BIT)
    && index_.is_valid()
    && index_ != Key::target_key()
    && !(is_arc() && !is_finite() );
}


inline void Combinatorial_curve::set_key(Key k) {
  index_=k;
}


inline void Combinatorial_curve::audit(unsigned int numv) const {
  if (is_special()) return;
  if (index_.is_input()) {
    CGAL_assertion(index_.input_index() < numv);
  } else {
    
  }
}


inline bool Combinatorial_curve::is_outward_rule() const {
  CGAL_precondition(is_rule());
  switch (rule_direction().is_outwards()) {
    return !is_inside();
  default:
    return is_inside();
  }
}
  
CGAL_AOS3_END_INTERNAL_NAMESPACE
