#include <CGAL/Constrained_triangulation_face_base_2.h>

namespace CGAL {

template <class Gt>
class Mesh_face_base_2 : public Constrained_triangulation_face_base_2<Gt>
{
public:
  typedef Gt Geom_traits;
  typedef Constrained_triangulation_face_base_2<Gt> CTFb;

protected:
  bool marked;

public:
  Mesh_face_base_2(): CTFb(), marked(false) {};

  Mesh_face_base_2(void* v0, void* v1, void* v2)
    : CTFb(v0,v1,v2), marked(false) {};

  Mesh_face_base_2(void* v0, void* v1, void* v2,
		   void* n0, void* n1, void* n2)
    : CTFb(v0,v1,v2,n0,n1,n2), marked(false) {};

  Mesh_face_base_2(void* v0, void* v1, void* v2,
		   void* n0, void* n1, void* n2,
		   bool c0, bool c1, bool c2 )
    : CTFb(v0,v1,v2,n0,n1,n2,c0,c1,c2), marked(false) {};

  inline
  bool is_marked() const { return marked; };

  inline
  void set_marked(const bool b) { marked=b; };
};

}; // namespace CGAL
