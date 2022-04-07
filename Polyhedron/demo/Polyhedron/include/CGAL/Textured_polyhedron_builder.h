#ifndef _TEXTURED_POLYHEDRON_BUILDER_
#define _TEXTURED_POLYHEDRON_BUILDER_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class Vertex_handle>
struct Vertex_cmp
{
  bool operator()(Vertex_handle a, Vertex_handle b) const
  {
    return &*a < &*b;
  }
};

template <class HDS,class Polyhedron,class Textured_polyhedron,class Kernel>
class Modifier_textured_polyhedron : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename Kernel::Point_3      Point;
  typedef typename HDS::Vertex          Vertex;
  typedef typename HDS::Face_handle     Face_handle;
  typedef typename HDS::Halfedge_handle Halfedge_handle;

  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;

  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
  Polyhedron *m_pMesh;
  std::map<Vertex_handle,int,Vertex_cmp<Vertex_handle> > m_vertex_map;

public:

  // life cycle
  Modifier_textured_polyhedron(Polyhedron *pMesh)
  {
    CGAL_assertion(pMesh != NULL);
    m_pMesh = pMesh;
  }
  ~Modifier_textured_polyhedron() {}

  void operator()( HDS& hds)
  {
    builder B(hds,true);
    B.begin_surface(3,1,6);
      add_vertices(B);
      add_facets(B);
    B.end_surface();
  }

  // add vertices
  void add_vertices(builder &B)
  {
    int index = 0;
    typename Polyhedron::Vertex_iterator it;
    for(it = m_pMesh->vertices_begin();
      it != m_pMesh->vertices_end();
      it++)
    {
      typename Polyhedron::Vertex_handle v = it;
      m_vertex_map[v] = index++;
      B.add_vertex(v->point());
    }
  }

  // add facets
  void add_facets(builder &B)
  {
    typename Polyhedron::Facet_iterator it;
    for(it = m_pMesh->facets_begin();
      it != m_pMesh->facets_end();
      it++)
    {
      B.begin_facet();
      HF_circulator he = it->facet_begin();
      HF_circulator end = he;
      CGAL_For_all(he,end)
        B.add_vertex_to_facet(m_vertex_map[he->vertex()]);
      B.end_facet();
    }
  }

};

template <class Polyhedron, class Textured_polyhedron, class Kernel>
class Textured_polyhedron_builder
{
public:
  typedef typename Textured_polyhedron::HalfedgeDS HalfedgeDS;
  Textured_polyhedron_builder() {}
  ~Textured_polyhedron_builder() {}

public:
  void run(Polyhedron &input,
    Textured_polyhedron &output)
  {
    Modifier_textured_polyhedron<HalfedgeDS,Polyhedron,Textured_polyhedron,Kernel> copier(&input);
    output.delegate(copier);
  }
};

#endif // _TEXTURED_POLYHEDRON_BUILDER_
