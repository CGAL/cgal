#ifndef CGAL_ARRANGEMENT_RULE_DIRECTION_H
#define CGAL_ARRANGEMENT_RULE_DIRECTION_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Coordinate_index.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
struct Combinatorial_curve;

struct Rule_direction {
  typedef Rule_direction This;
  typedef CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;
  explicit Rule_direction(int d);
  Rule_direction(): dir_(-1){};
  bool is_backwards() const;
  bool is_positive() const;
  bool is_negative() const;
  Coordinate_index constant_coordinate() const;
  bool is_vertical() const;
  bool can_intersect(const Combinatorial_curve &o) const;
  bool is_outwards() const;
  static Rule_direction make_from_part(int pt);
  CGAL_GETNR(int, part, return dir_);
  CGAL_COMPARISONS1(dir_);
  //CGAL_GETNR(Rule_direction, opposite,  return Rule_direction((dir_+2)%4));

  const char *to_str() const ;
  std::ostream& write(std::ostream& out) const{
    out << to_str();
    return out;
  }

  Rule_direction opposite() const;

  static Rule_direction right();
  static Rule_direction top();
  static Rule_direction left();
  static Rule_direction bottom();
  int index() const ;
private:
  int dir_;
};

CGAL_OUTPUT(Rule_direction);



CGAL_AOS3_END_INTERNAL_NAMESPACE


/*#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Rule_direction_impl.h>
#endif*/
#endif
