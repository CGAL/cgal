#ifndef CGAL_ARRANGEMENT_SPHERES_3_FEATURE_H
#define CGAL_ARRANGEMENT_SPHERES_3_FEATURE_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

class Rule_direction;

struct Combinatorial_curve{
  typedef Combinatorial_curve This;
  typedef CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;

  /* for rules
     Inside means below a horizontal arc and to the left of a vertical one
  */
  enum PART_BITS {R_BIT=2, T_BIT=4, L_BIT=8, B_BIT=16, ARC_BIT=32, /*INF_BIT=32,*/ IN_BIT=1};
  enum LOCATION_BITS {lR_BIT=R_BIT, lT_BIT=T_BIT, lL_BIT=L_BIT, lB_BIT=B_BIT, 
		      lIN_BIT=1, lOUT_BIT=32};
  enum Part {INVALID=0, L_RULE=L_BIT, R_RULE=R_BIT, T_RULE=T_BIT, B_RULE=B_BIT,
	     LB_ARC=L_BIT|B_BIT|ARC_BIT, 
	     LT_ARC=L_BIT|T_BIT|ARC_BIT,
	     RT_ARC=R_BIT|T_BIT|ARC_BIT, 
	     RB_ARC=R_BIT|B_BIT|ARC_BIT,
	     SPECIAL= R_BIT | L_BIT | B_BIT | T_BIT
	     /*, 
	       L_INF=L_BIT|INF_BIT,
	       R_INF=R_BIT|INF_BIT,
	       T_INF=T_BIT|INF_BIT, 
	       B_INF=B_BIT|INF_BIT*/};

  typedef Sphere_key Key;

  static Combinatorial_curve make_special(Key i) {
    return Combinatorial_curve(i, SPECIAL);
  }

  bool is_special() const {
    return pt_==SPECIAL;
  }

  void audit(unsigned int numvert) const;

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
  bool is_valid() const ;

  const Combinatorial_curve other_side() const;
  /*bool is_special() const {
    return pt_== SPECIAL;
    }*/

  bool is_inside() const {
    return pt_&IN_BIT;
  }
  void set_is_inside(bool tf) {
    if (tf) pt_= pt_| IN_BIT;
    else pt_= pt_ & (~IN_BIT);
  }

  const Combinatorial_curve strip_inside() const ;

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
  const Key& key() const {
    return index_;
  }
  void set_key(Key k);


  Coordinate_index
  constant_coordinate() const ;

  bool is_same_part(const Combinatorial_curve &o) const;

  bool can_intersect(int t) const;

  int quadrant() const;

  const Part part() const {
    //bool dont_use_part;
    return static_cast<Part>(pt_);
  }

  void flip_rule(Key k);

  bool is_rule() const {
    return ! is_arc();
  }
  bool is_negative() const {
    return pt_&L_BIT || pt_&B_BIT;
  }
  bool is_finite() const {
    return index_.is_input(); //! (pt_ & INF_BIT);
  }
  bool is_arc() const ;
  bool is_vertical() const ;

  bool is_same_side(Combinatorial_curve o) const ;

  static Rule_direction rule_direction(const Combinatorial_curve &a,
				       const Combinatorial_curve &b) ;

  CGAL_COMPARISONS2(index_, pt_);

  std::ostream &write(std::ostream&out) const;

  bool is_compatible_location(int i) const;



  Coordinate_index is_weakly_incompatible(int i) const ;

  Rule_direction rule_direction() const;
  int arc_index() const;

  static Combinatorial_curve make_rule(Key k, Rule_direction ruleindex);
  //static Combinatorial_curve make_rule(int ruleindex);

  static const char *to_string(int pt);

  const char *to_string() const {
    return to_string(pt_);
  }

  bool is_outward_rule() const ;
  
private:

  Combinatorial_curve(Key i, int pt, bool): index_(i), pt_(pt){
  }

  Key index_;
  int pt_;
};

CGAL_OUTPUT(Combinatorial_curve);

CGAL_AOS3_END_INTERNAL_NAMESPACE


/*#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve_impl.h>
#endif*/

#endif
