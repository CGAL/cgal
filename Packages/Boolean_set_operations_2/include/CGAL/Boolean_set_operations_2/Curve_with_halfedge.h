
#ifndef CURVE_WITH_HALFEDGE_H
#define CURVE_WITH_HALFEDGE_H

CGAL_BEGIN_NAMESPACE

template <class Arrangement_>
class Curve_with_halfedge
{
protected:
  typedef typename Arrangement_::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement_::Halfedge_const_handle  Halfedge_const_handle;
  
public:
  Halfedge_handle  m_he;

  Curve_with_halfedge()
  {};

  Curve_with_halfedge(Halfedge_handle he) : m_he(he)
  {}

  Halfedge_handle halfedge() const
  {
    return (m_he);
  }

  Halfedge_handle halfedge()
  {
    return (m_he);
  }

  void set_halfedge(Halfedge_handle he)
  {
    m_he = he;
  }
};

CGAL_END_NAMESPACE
#endif
