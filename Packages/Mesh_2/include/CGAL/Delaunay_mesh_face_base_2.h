#ifndef CGAL_DELAUNAY_FACE_BASE_2_H
#define CGAL_DELAUNAY_FACE_BASE_2_H

#include <CGAL/Constrained_triangulation_face_base_2.h>

namespace CGAL {

template <class Gt,
          class Fb = Constrained_triangulation_face_base_2<Gt> >
class Delaunay_mesh_face_base_2 : public Fb
{
public:
  typedef Gt Geom_traits;
  typedef Constrained_triangulation_face_base_2<Gt> CTFb;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename CTFb::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Delaunay_mesh_face_base_2<Gt,Fb2> Other;
  };

protected:
  bool marked;

public:
  Delaunay_mesh_face_base_2(): Fb(), marked(false) {};

  Delaunay_mesh_face_base_2(void* v0, void* v1, void* v2)
    : Fb(v0,v1,v2), marked(false) {};

  Delaunay_mesh_face_base_2(void* v0, void* v1, void* v2,
			    void* n0, void* n1, void* n2)
    : Fb(v0,v1,v2,n0,n1,n2), marked(false) {};

  Delaunay_mesh_face_base_2(void* v0, void* v1, void* v2,
			    void* n0, void* n1, void* n2,
			    bool c0, bool c1, bool c2 )
    : Fb(v0,v1,v2,n0,n1,n2,c0,c1,c2), marked(false) {};

  inline
  bool is_marked() const { return marked; };

  inline
  void set_marked(const bool b) { marked=b; };
};

}; // namespace CGAL

#endif
