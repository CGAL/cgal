#ifndef ARRANGEMENT_SPHERE_3_VERTEX_H
#define ARRANGEMENT_SPHERE_3_VERTEX_H

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>
#include <CGAL/basic.h>

/*
  If the two curves are both arcs, then the vertex is the one that is
  ccw with respect to the ordered centers.

  If they are both rules then the vertical one is first.
  If there is an arc and a rule then the rule is first.
*/
class Combinatorial_vertex{
public:
  typedef Combinatorial_curve::Key Key;
  enum Type_bits {R_BIT=1, SMALLER_BIT=2, SS_BIT=4, RR_BIT=8, SR_BIT=16};
  //enum Type {SS, SR, RR, SPECIAL};
  enum Type {SS=SS_BIT, RR= RR_BIT, RS=SR_BIT, 
	     SR=SR_BIT | R_BIT, INVALID=0,
	     SPECIAL=SS_BIT | RR_BIT | SR_BIT,
	     MAX = SR_BIT | RR_BIT|SS_BIT | SMALLER_BIT | R_BIT};

  Combinatorial_vertex(): type_(INVALID){}

  static Combinatorial_vertex make_special(Combinatorial_curve::Key i) ;



  static Combinatorial_vertex make_extremum(Key k,
					    Rule_direction dir) ;

  Combinatorial_vertex(Combinatorial_curve a, Combinatorial_curve b);

  void audit(unsigned int numvert) const;

  bool is_valid() const {
    //if (is_special()) return false;
    if (type_== INVALID) return false;
    CGAL_assertion(k_[0] != Key() && k_[1] != Key());
    return true;
  }

  bool operator==(const Combinatorial_vertex &o) const {
    return k_[0] == o.k_[0] && k_[1] == o.k_[1] && type_ == o.type_;
  }
  bool operator!=(const Combinatorial_vertex &o) const {
    return k_[0] != o.k_[0] || k_[1] != o.k_[1] || type_ == o.type_;
  }
  bool operator<(const Combinatorial_vertex &o) const {
    if (k_[0]<o.k_[0]) return true;
    else if (k_[0] > o.k_[0]) return false;
    else if (k_[1] < o.k_[1]) return true;
    else if (k_[1] > o.k_[1]) return false;
    else return (type_ < o.type_);
  }

  std::ostream &write(std::ostream &out) const;
  
  /*bool is_special() const {
    return a_.is_special() && b_.is_special();
    }*/
  /*Type type() const {
    if (b_.is_rule()) return RR;
    else if (a_.is_rule()) return SR;
    else return SS;
    }*/

  bool is_smaller() const {
    CGAL_precondition(is_sphere_rule());
    return type_&SMALLER_BIT;
  }

  bool is_sphere_sphere() const {
    return type_&SS_BIT;
  }

  bool is_sphere_rule() const {
    return type_&SR_BIT;
  }

  bool is_rule_rule() const {
    return type_& RR_BIT;
  }

  bool is_sphere_extremum() const {
    return is_sphere_rule() && k_[0] == k_[1];
  }

  bool is_special() const;
  
  Combinatorial_curve::Coordinate_index rule_coordinate() const {
    CGAL_precondition(is_sphere_rule());
    return plane_coordinate(type_&R_BIT);
  }

  Combinatorial_curve::Key sphere_key() const {
    CGAL_precondition(is_sphere_rule());
    if (type_&R_BIT) return k_[0];
    else return k_[1];
  }

  Combinatorial_curve::Key rule_key() const {
    CGAL_precondition(is_sphere_rule());
    if (type_&R_BIT) return k_[1];
    else return k_[0];
  }

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

  void set_key(Combinatorial_curve::Key k) {
    CGAL_precondition(is_sphere_extremum());
    k_[0]=k;
    k_[1]=k;
  }

  int extremum_index() const ;

  /*const Combinatorial_curve sphere(int i) const {
    CGAL_precondition(type()== SR || type()== SS);
    if (i==0 && type() == SS) {
    return a_;
    } else {
    return b_;
    }
    }

    const Combinatorial_curve rule(int i) const {
    if (i==0) {
    return a_;
    } else {
    CGAL_precondition( type()== RR);
    return b_;
    }
    }
  */

  /*void replace_rule(Combinatorial_curve c) {
    if (is_sphere_rule()) {
      CGAL_assertion(0);
      //CGAL_precondition(c.is_vertical() == a_.is_vertical());
      CGAL_precondition(!is_sphere_extremum());
      k_[c.constant_coordinate()]= c.key();
    } else if (is_rule_rule()) {
      if (c.is_vertical() && a_.is_vertical()){
	a_=c;
      } else {
	b_=c;
      }
    } else {
      CGAL_assertion(0);
    }
    }*/

  /*const Combinatorial_curve other_curve(Combinatorial_curve c) const {
    if (is_rule_rule()) {
      CGAL_precondition(c.is_rule());
      bool do_i
      if (c.constant_coordinate()==0) return b_;
      if (c.constant_coordinate()==1) return a_;
    } else if (is_sphere_rule()) {
      if (c.is_rule()) return sphere(0);
      else return rule(0);
    } else {
      if (c== sphere(0)) return sphere(1);
      else return sphere(0);
    }
    }*/

  

  /*const Combinatorial_curve first() const {
    return a_;
  }

  const Combinatorial_curve second() const {
    return b_;
  }
  Combinatorial_curve& first()  {
    return a_;
  }
  Combinatorial_curve& second()  {
    return b_;
    }*/
  

  bool is_finite() const {
    return k_[0].is_input() && k_[1].is_input();
  }

protected:
  Combinatorial_curve::Key k_[2];
  Type type_;
};

inline std::ostream &operator<<(std::ostream &out, Combinatorial_vertex f) {
  return f.write(out);
}

#endif
