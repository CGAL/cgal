
#ifndef POINT_WITH_VERTEX_H
#define POINT_WITH_VERTEX_H

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Point_with_vertex
{
protected:
  typedef typename Arrangement_::Vertex_handle        Vertex_handle;
  typedef typename Arrangement_::Vertex_const_handle  Vertex_const_handle;
  
public:
  Vertex_handle  m_v;

  Point_with_vertex()
  {};

  Point_with_vertex(Vertex_handle v) : m_v(v)
  {}

  Vertex_handle vertex() const
  {
    return (m_v);
  }

  Vertex_handle vertex()
  {
    return (m_v);
  }

  void set_halfedge(Vertex_handle v)
  {
    m_v = v;
  }
};

CGAL_END_NAMESPACE
#endif
