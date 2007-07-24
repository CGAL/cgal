#ifndef ARRANGEMENT_SPHERE_3_VERTEX_H
#define ARRANGEMENT_SPHERE_3_VERTEX_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3/coordinates.h>
#include <CGAL/basic.h>
CGAL_BEGIN_NAMESPACE


namespace Arrangement_of_spheres_3_internal {

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

    CGAL_KINETIC_IS(smaller,  
		    CGAL_precondition(is_sphere_rule());
		    return type_&SMALLER_BIT);
    CGAL_KINETIC_IS(sphere_sphere,
		    return type_&SS_BIT);
    CGAL_KINETIC_IS(sphere_rule,
		    return type_&SR_BIT);
    CGAL_KINETIC_IS(rule_rule,
		    return type_& RR_BIT);
    CGAL_KINETIC_IS(sphere_extremum,
		    return is_sphere_rule() && k_[0] == k_[1]);

    bool is_special() const;
  
    CGAL_KINETIC_GETNR(Combinatorial_curve::Coordinate_index,
			    rule_coordinate,
			    CGAL_precondition(is_sphere_rule());
			    return plane_coordinate(type_&R_BIT));

    CGAL_KINETIC_GET(Combinatorial_curve::Key,
			  sphere_key,
			  CGAL_precondition(is_sphere_rule());
			  if (type_&R_BIT) return k_[0];
			  else return k_[1]);

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


    void set_rule_key(Combinatorial_curve::Key k) {
      CGAL_precondition(is_sphere_rule());
      if (type_&R_BIT) k_[1]=k;
      else k_[0]=k;
    }

    void swap_key(Combinatorial_curve::Key a, Combinatorial_curve::Key b) {
      if (k_[0]==a) k_[0]=b;
      if (k_[1]==a) k_[1]=b;
    }

    Combinatorial_curve::Key other_key(Combinatorial_curve::Key a) const {
      if (k_[0]== a) return k_[1];
      else return k_[0];
    }

    Coordinate_index other_curve_constant_coordinate(Combinatorial_curve a) const {
      CGAL_precondition(is_sphere_rule() || is_rule_rule());
      if (a.is_arc()) {
	return rule_coordinate();
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

  protected:
    Combinatorial_curve::Key k_[2];
    Type type_;
  };

  inline std::ostream &operator<<(std::ostream &out, Combinatorial_vertex f) {
    return f.write(out);
  }

}
CGAL_END_NAMESPACE

#ifdef CGAL_ARRANGEMENT_OF_SPHERES_3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex_impl.h>
#endif

#endif
