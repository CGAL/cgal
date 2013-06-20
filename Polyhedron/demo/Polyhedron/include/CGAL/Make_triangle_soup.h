#ifndef _MAKE_SOUP_
#define _MAKE_SOUP_

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

template<class Datum, class Polyhedron, class Kernel>
struct CModifierTriangleSoupVertexMapper 
{ };

template<class Polyhedron, class Kernel>
struct CModifierTriangleSoupVertexMapper<typename Kernel::Triangle_3, Polyhedron, Kernel>
{
  typedef typename Kernel::Triangle_3 Datum;

  const typename Kernel::Point_3& operator()(const Datum& datum, unsigned int index) const {
    return datum[index];
  }
};

template<class Polyhedron, class Kernel>
struct CModifierTriangleSoupVertexMapper<typename Polyhedron::Facet_const_handle, Polyhedron, Kernel>
{
  typedef typename Polyhedron::Facet_const_handle Datum;

  const typename Kernel::Point_3& operator()(const Datum& datum, unsigned int index) const{
    typename Polyhedron::Halfedge_around_facet_const_circulator circ= datum->facet_begin();
    std::advance(circ, index);
    return circ->vertex()->point();
  }
};

template <class HDS,class Polyhedron,class Kernel,class InputIterator>
class CModifierTriangleSoup : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
  InputIterator m_begin, m_end;

public:

  // life cycle
  CModifierTriangleSoup(InputIterator begin,InputIterator end)
    : m_begin(begin), m_end(end)
  {
  }
  ~CModifierTriangleSoup() {}

  // make a polygon soup
  void operator()( HDS& hds)
  {
    CModifierTriangleSoupVertexMapper<
     typename std::iterator_traits<InputIterator>::value_type, Polyhedron, Kernel> mapper;

    builder B(hds,true);
    B.begin_surface(3,1,6);

    int index = 0;
    InputIterator it;
    for(it = m_begin; it != m_end; it++)
    {
      B.add_vertex(mapper(*it, 0));
      B.add_vertex(mapper(*it, 1));
      B.add_vertex(mapper(*it, 2));

      B.begin_facet();
      B.add_vertex_to_facet(index++);
      B.add_vertex_to_facet(index++);
      B.add_vertex_to_facet(index++);
      B.end_facet();
    }
    B.end_surface();
  }
};

template <class Polyhedron,class Kernel,class InputIterator>
class Make_triangle_soup
{
public:
  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  Make_triangle_soup() {}
  ~Make_triangle_soup() {}

public:
  void run(InputIterator begin,
    InputIterator end,
    Polyhedron &output)
  {
    CModifierTriangleSoup<HalfedgeDS,Polyhedron,Kernel,InputIterator> soup(begin,end);
    output.delegate(soup);
  }
};

#endif // _MAKE_SOUP_
