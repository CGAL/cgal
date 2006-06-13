#ifndef CGAL_SPHERE_KEY_H
#define CGAL_SPHERE_KEY_H
struct Sphere_key{
  static const int BL=-2, TR=-1, TEMP=-3;
  Sphere_key(): id_(-4){}
  explicit Sphere_key(int i): id_(i){}
  int to_index() const {return id_;}
  bool is_valid() const {
    return id_ > -4;
  }
  bool operator==(const Sphere_key &o) const {
    return id_==o.id_;
  }
  bool operator!=(const Sphere_key &o) const {
    return id_!=o.id_;
  }
  bool operator<(const Sphere_key &o) const {
    return id_<o.id_;
  }
  bool operator>(const Sphere_key &o) const {
    return id_>o.id_;
  }
  bool is_input() const {
    return id_>=0;
  }
  unsigned int input_index() const {
    CGAL_precondition(is_input());
    return id_;
  }

  unsigned int internal_index() const {
    return id_+3;
  }

  std::ostream &write(std::ostream &out) const {
    if (id_==TEMP) out << "<T>";
    else if (id_ == TR) out << "<TR>";
    else if (id_ == BL) out << "<BL>";
    else out << "<" << id_ << ">";
    return out;
  }

  static Sphere_key temp_key() {return Sphere_key(TEMP);}
  static Sphere_key bl_key() {return Sphere_key(BL);}
  static Sphere_key tr_key() {return Sphere_key(TR);}

  int id_;
};

inline std::ostream &operator<<(std::ostream &out, Sphere_key sk) {
  return sk.write(out);
}
#endif
