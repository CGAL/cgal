#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>

Combinatorial_vertex::Combinatorial_vertex(Combinatorial_curve a, 
					   Combinatorial_curve b){
  CGAL_precondition(a.key() != b.key() 
		    || !a.key().is_input() && !b.key().is_input());
  if (a.is_rule() && b.is_rule()){
    k_[project(a.constant_coordinate())]= a.key();
    k_[project(b.constant_coordinate())]= b.key();
    type_=RR;
  } else if (a.is_rule()){
    k_[project(a.constant_coordinate())]=a.key();
    k_[project(other_plane_coordinate(a.constant_coordinate()))]= b.key();
    type_= static_cast<Type>(SR_BIT + project(a.constant_coordinate()));
    if (!a.is_negative() && !a.is_same_side(b)
	|| a.is_negative() && a.is_same_side(b)) {
      type_ = static_cast<Type>(type_ | SMALLER_BIT);
    }
  } else if (b.is_rule()) {
    k_[project(b.constant_coordinate())]=b.key();
    k_[project(other_plane_coordinate(b.constant_coordinate()))]= a.key();
    type_= static_cast<Type>(SR_BIT + project(b.constant_coordinate()));
    if (!b.is_negative() && !a.is_same_side(b)
	|| b.is_negative() && a.is_same_side(b)) {
      type_ = static_cast<Type>(type_ | SMALLER_BIT);
    }
  } else {
    k_[0]=a.key();
    k_[1]=b.key();
    type_= SS;
  }	
  CGAL_assertion(k_[0] != k_[1]
		 || !k_[0].is_input() && !k_[1].is_input());
  CGAL_postcondition(is_valid());
}

Combinatorial_vertex 
Combinatorial_vertex::make_special(Combinatorial_curve::Key i) {
  Combinatorial_vertex cv;
  cv.type_= SPECIAL;
  cv.k_[0]=i;
  cv.k_[1]=i;
  return cv;
}


bool Combinatorial_vertex::is_special() const {
  return type_== SPECIAL || k_[0].is_target();
}


Combinatorial_vertex Combinatorial_vertex::make_extremum(Key k,
								int dir) {
  Combinatorial_vertex ret;
  ret.k_[0]=k;
  ret.k_[1]=k;
  if (dir%2==0) {
    ret.type_= SR;
  } else {
    ret.type_= RS;
  }
  if (dir >= 2) {
    ret.type_= static_cast<Type>(ret.type_| SMALLER_BIT);
  }
  return ret;
}

int Combinatorial_vertex::extremum_index() const {
  CGAL_precondition(is_sphere_extremum());
  if (rule_coordinate() == plane_coordinate(1)){
    if (!is_smaller()) return 0;
    else return 2;
  } else {
    if (!is_smaller()) return 1;
    else return 3;
  }
}

std::ostream &Combinatorial_vertex::write(std::ostream &out) const {
  if ( type_ & RR_BIT) {
    out << k_[0] << "," << k_[1];
  } else if (is_sphere_extremum()) {
    out << k_[0];
    switch(extremum_index()) {
    case 0:
      out << "R"; break;
    case 1:
      out << "T"; break;
    case 2: 
      out << "L"; break;
    default:
      out << "B"; break;
    }
  } else if (type_ & SR_BIT) {
   
    if (rule_coordinate() == plane_coordinate(0)){
      out << rule_key() << "," << sphere_key() << "s";
    } else {
      out << sphere_key() << "s," << rule_key();
    }
    
  } else if (type_ &SS_BIT) {
    out << k_[0] << ":" << k_[1];
  } else if (type_ == SPECIAL) {
    out << "SPECIAL";
  } else {
    out << "INVALID";
  }
  return out;
}


bool Combinatorial_curve::can_intersect(int t) const {
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


int Combinatorial_curve::quadrant() const {
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


int Combinatorial_curve::rule_direction(const Combinatorial_curve &a,
					       const Combinatorial_curve &b) {
  CGAL_precondition(a.is_arc());
  CGAL_precondition(b.is_arc());
  CGAL_precondition(a.key() == b.key());
  switch(a.pt_ & b.pt_ & (L_BIT | T_BIT | R_BIT | B_BIT)) {
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
  }
}
  

bool Combinatorial_curve::is_compatible_location(int i) const {
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

const char *Combinatorial_curve::to_string(int pt){
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


bool Combinatorial_curve::is_arc() const {
  bool nb= (pt_ &ARC_BIT) != 0;
  if (pt_&ARC_BIT) {
    CGAL_assertion(nb);
  }
  if (!(pt_&ARC_BIT)) {
    CGAL_assertion(!nb);
  }
  return nb;
}
bool Combinatorial_curve::is_vertical() const {
  CGAL_precondition(is_rule());
  //if (is_finite()) {
  return (pt_ & T_BIT) || (pt_ & B_BIT);
  /*} else {
    return (pt_ & L_BIT) || (pt_ & R_BIT);
    }*/
}
bool Combinatorial_curve::is_same_side(Combinatorial_curve o) const {
  int u= o.pt_ & pt_; 
  return u &( T_BIT | L_BIT | R_BIT | B_BIT);
}


std::ostream &Combinatorial_curve::write(std::ostream&out) const {
  CGAL_precondition(is_valid());
  out << to_string(pt_);
  /*if (is_finite())*/  out << index_;
  return out;
}

bool Combinatorial_curve::can_intersect(const Combinatorial_curve &o) const {
  CGAL_assertion(is_rule());
  return pt_&o.quadrant();
}

Coordinate_index Combinatorial_curve::is_weakly_incompatible(int i) const {
  int a= i&pt_&(~lOUT_BIT) &(~lIN_BIT);
  if (a== L_BIT || a==R_BIT) return plane_coordinate(1);
  else if (a== T_BIT || a== B_BIT) return plane_coordinate(0);
  else return Coordinate_index();
}

int Combinatorial_curve::rule_index() const {
  CGAL_precondition(is_rule());
  if (pt_&R_RULE) return 0;
  else if (pt_&T_RULE) return 1;
  else if (pt_&L_RULE) return 2;
  else return 3;
}

int Combinatorial_curve::arc_index() const {
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

Combinatorial_curve Combinatorial_curve::make_rule(Key k, int ruleindex) {
  switch(ruleindex) {
  case 0:
    return Combinatorial_curve(k, R_RULE);
  case 1:
    return Combinatorial_curve(k, T_RULE);
  case 2:
    return Combinatorial_curve(k, L_RULE);
  default:
    return Combinatorial_curve(k, B_RULE);
  }
}


Coordinate_index
Combinatorial_curve::constant_coordinate() const {
  CGAL_precondition(is_rule());
  Coordinate_index r;
  if ((pt_&R_BIT) || (pt_ &L_BIT)) 
    r= plane_coordinate(1);
  else r= plane_coordinate(0);
  //if (!is_finite()) r=1-r;
  return r;
}

bool Combinatorial_curve::is_same_part(const Combinatorial_curve &o) const {
  CGAL_precondition(is_inside() == o.is_inside());
  CGAL_precondition(is_rule());
  CGAL_precondition(o.is_rule());
  return o.pt_==pt_;
}

void Combinatorial_curve::flip_rule(Key k) {
  CGAL_precondition(is_rule());
  index_=k;
  if (is_vertical()) {
    pt_ = pt_^ (T_BIT | B_BIT);
  } else {
    pt_ = pt_^ (L_BIT | R_BIT);
  }
}


const Combinatorial_curve Combinatorial_curve::strip_inside() const {
  Combinatorial_curve ret(index_, pt_& (~IN_BIT), true);
  CGAL_assertion(!ret.is_inside());
  CGAL_assertion(ret.is_valid());
  return ret;
}

const Combinatorial_curve Combinatorial_curve::other_side() const {
  Combinatorial_curve ret=*this;
  ret.set_is_inside(!is_inside());
  CGAL_assertion(!(is_inside() && ret.is_inside()));
  return ret;
}


bool Combinatorial_curve::is_valid() const {
  return pt_>0 && pt_ < 128
    //&& (pt_ ^ IN_BIT)
    && index_.is_valid()
    && index_ != Key::target_key()
    && !(is_arc() && !is_finite() );
}


void Combinatorial_curve::set_key(Key k) {
  index_=k;
}
