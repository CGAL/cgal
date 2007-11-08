#ifndef ARRANGEMENT_SPHERE_3_VERTEX_H
#define ARRANGEMENT_SPHERE_3_VERTEX_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/basic.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


/*
  If the two curves are both arcs, then the vertex is the one that is
  ccw with respect to the ordered centers.

  If they are both rules then the vertical one is first.
  If there is an arc and a rule then the rule is first.
*/
class Combinatorial_vertex{
  typedef Combinatorial_vertex This;
public:
  typedef Combinatorial_curve::Key Key;
  enum Type_bits {R_BIT=1, SMALLER_BIT=2, SS_BIT=4, RR_BIT=8, SR_BIT=16, RULE_SIGN=32};
  //enum Type {SS, SR, RR, SPECIAL};
  enum Type {SS=SS_BIT, RR= RR_BIT, RS=SR_BIT, 
	     SR=SR_BIT | R_BIT, INVALID=0,
	     SPECIAL=SS_BIT | RR_BIT | SR_BIT,
	     MAX = SR_BIT | RR_BIT|SS_BIT | SMALLER_BIT | R_BIT};

  Combinatorial_vertex(): type_(INVALID){}


  static Combinatorial_vertex make_extremum(Key k,
					    Rule_direction dir) ;


 
  Combinatorial_vertex(Combinatorial_curve a, Combinatorial_curve b);

  static Combinatorial_vertex make_rule_rule(Key a, Key b) {
    Combinatorial_vertex ret;
    ret.k_[0]=a;
    ret.k_[1]=b;
    ret.type_= RR;
    return ret;
  }

  void audit(unsigned int numvert) const;

  bool is_valid() const {
    //if (is_special()) return false;
    if (type_== INVALID) return false;
    CGAL_assertion(k_[0] != Key() && k_[1] != Key());
    return true;
  }

  CGAL_COMPARISONS3(k_[0], k_[1], type_);


  std::ostream &write(std::ostream &out) const;
  
  /*bool is_special() const {
    return a_.is_special() && b_.is_special();
    }*/
  /*Type type() const {
    if (b_.is_rule()) return RR;
    else if (a_.is_rule()) return SR;
    else return SS;
    }*/

 
  CGAL_IS(sphere_sphere,
	  return type_&SS_BIT);
  CGAL_IS(sphere_rule,
	  return type_&SR_BIT);
  CGAL_IS(rule_rule,
	  return type_& RR_BIT);
  CGAL_IS(sphere_extremum,
	  return is_sphere_rule() && k_[0] == k_[1]);

  bool is_special() const;
  
  CGAL_GETNR(Combinatorial_curve::Coordinate_index,
	     rule_constant_coordinate,
	     CGAL_precondition(is_sphere_rule());
	     return plane_coordinate(type_&R_BIT));
  
  CGAL_GET(Combinatorial_curve::Key,
	   sphere_key,
	   CGAL_precondition(is_sphere_rule());
	   if (type_&R_BIT) return k_[0];
	   else return k_[1]);
  
  Combinatorial_curve::Key rule_key() const {
    CGAL_precondition(is_sphere_rule());
    if (type_&R_BIT) return k_[1];
    else return k_[0];
  }

  /*bool rule_is_positive() const {
    CGAL_precondition(is_sphere_rule());
    return !(type__&SMALLER_BIT);
    }*/

  /*Rule_direction rule_direction() const {
    CGAL_precondition(is_sphere_rule());
    if (is_smaller() && rule_constant_coordinate() == plane_coordinate(0)) {
    return Rule_direction(0);
    } else if (!is_smaller() 
    && rule_constant_coordinate() == plane_coordinate(0)) {
    return Rule_direction(2);
    } else if (is_smaller() 
    && rule_constant_coordinate() == plane_coordinate(1)) {
    return Rule_direction(1);
    } else if (!is_smaller() 
    && rule_constant_coordinate() == plane_coordinate(1)) {
    return Rule_direction(3);
    } else {
    CGAL_error();
    return Rule_direction();
    }
    }*/

  Combinatorial_curve::Key key() const {
    CGAL_precondition(is_sphere_extremum());
    return k_[0];
  }

  Combinatorial_curve::Key rule_key(Coordinate_index i) const {
    CGAL_precondition(is_rule_rule());
    CGAL_precondition(i == plane_coordinate(0) || i == plane_coordinate(1)); 
    return k_[project(i)];
  }

  Combinatorial_curve::Key sphere_key(int i) const {
    CGAL_precondition(is_sphere_sphere());
    return k_[i];
  }

  void set_rule_key(Combinatorial_curve::Key k) {
    CGAL_precondition(is_sphere_rule());
    if (type_&R_BIT) k_[1]=k;
    else k_[0]=k;
  }

 
  Combinatorial_curve::Key other_key(Combinatorial_curve::Key a) const {
    if (k_[0]== a) return k_[1];
    else return k_[0];
  }


  bool other_curve_key(Combinatorial_curve a) const {
    if (a.is_arc()) {
      CGAL_precondition(is_sphere_rule() || is_sphere_sphere());
      return is_sphere_rule();
    } else {
      CGAL_precondition(is_sphere_rule() || is_rule_rule());
      return is_rule_rule();
    }
  }

  

  Rule_direction sphere_extremum_index() const ;

  

  bool is_finite() const {
    return k_[0].is_input() && k_[1].is_input();
  }

  // Return true if the line defining the rule which defines this point
  // points in the positive direction (and so it is the smaller of the two points)
  CGAL_IS(smaller,  
	  CGAL_precondition(is_sphere_rule());
	  return type_&SMALLER_BIT);

  
  void swap_key(Combinatorial_curve::Key a, Combinatorial_curve::Key b) {
    if (k_[0]==a) k_[0]=b;
    if (k_[1]==a) k_[1]=b;
  }

protected:

  void set_key(Combinatorial_curve::Key k) {
    CGAL_precondition(is_sphere_extremum());
    k_[0]=k;
    k_[1]=k;
  }



  Coordinate_index other_curve_constant_coordinate(Combinatorial_curve a) const {
    CGAL_precondition(is_sphere_rule() || is_rule_rule());
    if (a.is_arc()) {
      return rule_constant_coordinate();
    } else {
      CGAL_precondition(is_rule_rule());
      return other_plane_coordinate(a.constant_coordinate());
    }
  }

  bool other_curve_is_rule(Combinatorial_curve a) const {
    if (a.is_arc()) {
      CGAL_precondition(is_sphere_rule() || is_sphere_sphere());
      return is_sphere_rule();
    } else {
      CGAL_precondition(is_sphere_rule() || is_rule_rule());
      return is_rule_rule();
    }
  }


  static Combinatorial_vertex make_special(Combinatorial_curve::Key i) ;


 
  Combinatorial_curve::Key k_[2];
  Type type_;
};


CGAL_OUTPUT(Combinatorial_vertex);
  
CGAL_AOS3_END_INTERNAL_NAMESPACE


/*#ifdef CGAL_AOS3_USE_TEMPLATES
  #include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex_impl.h>
  #endif*/

#endif
