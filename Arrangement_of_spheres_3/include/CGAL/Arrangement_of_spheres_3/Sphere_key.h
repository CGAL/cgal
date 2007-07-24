#ifndef CGAL_SPHERE_KEY_H
#define CGAL_SPHERE_KEY_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

struct Sphere_key{
  typedef Sphere_key This;
  enum Labels {BL=-2, TR=-1, TEMP=-3, TARGET=-5, INVALID=-4};
  Sphere_key(): id_(INVALID){}
  explicit Sphere_key(int i): id_(i){}
  CGAL_GETNR(int, index,return id_);
  CGAL_GETNR(unsigned int, input_index,CGAL_precondition(is_input()); return id_);
  CGAL_GETNR(unsigned int, internal_index,return id_+3);
  bool is_valid() const {
    return id_ > -4;
  }
  CGAL_COMPARISONS1(id_);
  CGAL_IS(input, return id_>=0);
  CGAL_IS(target, return id_==TARGET);
 

  std::ostream &write(std::ostream &out) const {
    if (id_==TARGET) out << "tar";
    else if (id_ == TR) out << "tr";
    else if (id_ == BL) out << "bl";
    else out << id_;
    return out;
  }

  static Sphere_key target_key() {return Sphere_key(TARGET);}
  static Sphere_key temp_key() {return Sphere_key(TEMP);}
  static Sphere_key bl_key() {return Sphere_key(BL);}
  static Sphere_key tr_key() {return Sphere_key(TR);}

  int id_;
};

CGAL_OUTPUT(Sphere_key);
CGAL_AOS3_END_INTERNAL_NAMESPACE
#endif
