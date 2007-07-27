#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE



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
							 Rule_direction dir) {
  Combinatorial_vertex ret;
  ret.k_[0]=k;
  ret.k_[1]=k;
  if (!dir.is_vertical()) {
    ret.type_= SR;
  } else {
    ret.type_= RS;
  }
  if (dir.is_negative()) {
    ret.type_= static_cast<Type>(ret.type_| SMALLER_BIT);
  }
  return ret;
}

 Rule_direction Combinatorial_vertex::sphere_extremum_index() const {
  CGAL_precondition(is_sphere_extremum());
  if (rule_coordinate() == plane_coordinate(1)){
    if (!is_smaller()) return Rule_direction(0);
    else return Rule_direction(2);
  } else {
    if (!is_smaller()) return Rule_direction(1);
    else return Rule_direction(3);
  }
}

 void Combinatorial_vertex::audit(unsigned int numv) const {
  if (is_special()) return;
  if (k_[0].is_input()) {
    CGAL_assertion(k_[0].input_index() < numv);
  }
  if (k_[1].is_input()) {
    CGAL_assertion(k_[1].input_index() < numv);
  }
}

 std::ostream &Combinatorial_vertex::write(std::ostream &out) const {
  if ( type_ & RR_BIT) {
    out << k_[0] << "," << k_[1];
  } else if (is_sphere_extremum()) {
    out << k_[0];
    Rule_direction d= sphere_extremum_index();
    if (d== Rule_direction(0)) {
      out << "R";
    } else if (d== Rule_direction(1)) {
      out << "T";
    } else if (d== Rule_direction(2)) {
      out << "L";
    } else {
      out << "B";
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

CGAL_AOS3_END_INTERNAL_NAMESPACE
