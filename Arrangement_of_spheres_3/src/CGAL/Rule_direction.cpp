#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
#define inline 

inline Rule_direction::Rule_direction(int d): dir_(2){
  //CGAL_precondition(d >=0 && d < 4);
  dir_= dir_ << d;
  if (d==0) CGAL_assertion(dir_== Combinatorial_curve::R_BIT);
  else if (d==1) CGAL_assertion(dir_==  Combinatorial_curve::T_BIT);
  else if (d==2) CGAL_assertion(dir_== Combinatorial_curve::L_BIT);
  else if (d==3) CGAL_assertion(dir_==  Combinatorial_curve::B_BIT);
  else CGAL_assertion(0);
}


inline bool Rule_direction::is_backwards() const {
  return dir_== Combinatorial_curve::T_BIT 
    || dir_==Combinatorial_curve::L_BIT;
}



inline bool Rule_direction::is_positive() const {
  return dir_== Combinatorial_curve::T_BIT ||
    dir_==Combinatorial_curve::R_BIT;
}


inline bool Rule_direction::is_negative() const {
  return dir_== Combinatorial_curve::B_BIT ||
    dir_==Combinatorial_curve::L_BIT;
}


inline Coordinate_index Rule_direction::constant_coordinate() const {
  if (!is_vertical()) return plane_coordinate(1);
  else return plane_coordinate(0);
}


inline bool Rule_direction::is_vertical() const {
  return (dir_== Combinatorial_curve::T_BIT 
	  || dir_== Combinatorial_curve::B_BIT );
}


inline bool Rule_direction::can_intersect(const Combinatorial_curve &o) const{
  return dir_&o.quadrant();
}


inline bool Rule_direction::is_outwards() const {
  return dir_ == Combinatorial_curve::R_BIT
    || dir_== Combinatorial_curve::B_BIT;
}


inline const char *Rule_direction::to_str() const {
  switch( dir_) {
  case(Combinatorial_curve::T_BIT): return "T";
  case(Combinatorial_curve::B_BIT): return "B";
  case(Combinatorial_curve::L_BIT): return "L";
  case(Combinatorial_curve::R_BIT): return "R";
  default: return "INV";
  }
}


inline int Rule_direction::index() const {
  switch (dir_) {
  case(Combinatorial_curve::T_BIT): return 1;
  case(Combinatorial_curve::B_BIT): return 3;
  case(Combinatorial_curve::L_BIT): return 2;
  case(Combinatorial_curve::R_BIT): return 0;
  default: 
    CGAL_assertion(0);
    return -1;
  }
}



inline Rule_direction Rule_direction::right(){
  return make_from_part(Combinatorial_curve::R_BIT);
}



inline Rule_direction Rule_direction::top(){
  return make_from_part(Combinatorial_curve::T_BIT);
}



inline Rule_direction Rule_direction::left(){
  return make_from_part(Combinatorial_curve::L_BIT);
}


inline Rule_direction Rule_direction::bottom(){
  return make_from_part(Combinatorial_curve::B_BIT);
}



inline Rule_direction Rule_direction::make_from_part(int pt) {
  CGAL_assertion(pt <= ( Combinatorial_curve::B_BIT) && pt >0);
  Rule_direction r;
  r.dir_=pt;
  if (pt & Combinatorial_curve::IN_BIT) {
    return r.opposite();
  } else return r;
}

inline Rule_direction Rule_direction::opposite() const {
  switch( dir_) {
  case(Combinatorial_curve::T_BIT): 
    return make_from_part(Combinatorial_curve::B_BIT);
  case(Combinatorial_curve::B_BIT): 
    return make_from_part(Combinatorial_curve::T_BIT);
  case(Combinatorial_curve::L_BIT): 
    return make_from_part(Combinatorial_curve::R_BIT);
  case(Combinatorial_curve::R_BIT): 
    return make_from_part(Combinatorial_curve::L_BIT);
  default: 
    CGAL_assertion(0);
    return Rule_direction();
  }
}
CGAL_AOS3_END_INTERNAL_NAMESPACE
