#ifndef _MAKE_SOUP_
#define _MAKE_SOUP_

#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

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
    builder B(hds,true);
    B.begin_surface(3,1,6);

    int index = 0;
    InputIterator it;
    for(it = m_begin; it != m_end; it++)
    {
      const Triangle& triangle = *it;

      B.add_vertex(triangle[0]);
      B.add_vertex(triangle[1]);
      B.add_vertex(triangle[2]);

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
