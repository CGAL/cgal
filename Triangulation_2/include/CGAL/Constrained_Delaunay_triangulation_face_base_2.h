#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

namespace CGAL {

template <typename Kernel,
          class Fb = CGAL::Constrained_triangulation_face_base_2<Kernel> >
class Constrained_Delaunay_triangulation_face_base_2
  : public Fb
{
public:
  typedef Fb Base;
  typedef typename Base::Vertex_handle  Vertex_handle;
  typedef typename Base::Face_handle    Face_handle;
  typedef std::pair<Face_handle, int> Edge;
  typedef Constrained_Delaunay_triangulation_face_base_2 CDT_face_base;

  typedef typename Kernel::Point_2      Point;
  typedef typename Kernel::Segment_2    Segment;

private:
  bool m_blind;
  Edge m_blinding_constraint;

public:
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Base::template Rebind_TDS<TDS2>::Other Fb2;
    typedef Constrained_Delaunay_triangulation_face_base_2<Kernel,Fb2>
      Other;
  };

public:
  CDT_face_base()
    : Base(),
      m_blind(false)
  {
  }
  CDT_face_base(Vertex_handle v1,
                Vertex_handle v2,
                Vertex_handle v3)
    : Base(v1,v2,v3),
      m_blind(false)
  {
  }
  CDT_face_base(Vertex_handle v1,
                Vertex_handle v2,
                Vertex_handle v3,
                Face_handle f1,
                Face_handle f2,
                Face_handle f3)
    : Base(v1,v2,v3,f1,f2,f3),
      m_blind(false)
  {
  }
  CDT_face_base(Face_handle f)
    : Base(f),
      m_blind(false)
  {
  }

  // sees its circumcenter or not?
  const bool& blind() const { return m_blind; }
  bool& blind(){ return m_blind; }

  // if blind, the constrained edge that prevents the face
  // to see its circumcenter 
  const Edge& blinding_constraint() const { return m_blinding_constraint; }
  Edge& blinding_constraint() { return m_blinding_constraint; }
};

} //namespace CGAL 

#endif //CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_FACE_BASE_2_2_H
