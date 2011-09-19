#ifndef _DUALIZER_
#define _DUALIZER_

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <class Face_handle>
struct Facet_cmp
{
  bool operator()(Face_handle a, Face_handle b) const
  {
    return &*a < &*b;
  }
};

template <class HDS,class Polyhedron,class Kernel>
class CModifierDual : public CGAL::Modifier_base<HDS>
{
private:
  typedef typename Kernel::Point_3      Point;
  typedef typename Kernel::Plane_3      Plane;
  typedef typename Kernel::Vector_3     Vector;
  typedef typename Kernel::FT           FT;

  typedef typename HDS::Vertex          Vertex;
  typedef typename HDS::Face_handle     Face_handle;
  typedef typename HDS::Halfedge_handle Halfedge_handle;

  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator HV_circulator;

  typedef typename CGAL::Polyhedron_incremental_builder_3<HDS> builder;
  Polyhedron *m_pMesh;
  std::map<Face_handle,int,Facet_cmp<Face_handle> > m_face_map;

public:

  // life cycle
  CModifierDual(Polyhedron *pMesh)
  {
    CGAL_assertion(pMesh != NULL);
    m_pMesh = pMesh;
  }
  ~CModifierDual() {}

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
    Facet_iterator it;
    for(it = m_pMesh->facets_begin();
      it != m_pMesh->facets_end();
      it++)
    {
      Face_handle f = it;
      m_face_map[f] = index++;
      B.add_vertex(dual(f));
    }
  }

  Plane facet_plane(Face_handle f)
  {
    const Point& a = f->halfedge()->vertex()->point();
    const Point& b = f->halfedge()->next()->vertex()->point();
    const Point& c = f->halfedge()->next()->next()->vertex()->point();
    return Plane(a,b,c);
  }

  Point dual(Face_handle f)
  {
    Plane plane = facet_plane(f);
    FT sqd = CGAL::squared_distance(Point(CGAL::ORIGIN),plane);
    FT distance_to_origin = std::sqrt(sqd);
    Vector normal = plane.orthogonal_vector();
    normal = normal / std::sqrt(normal * normal);
    return CGAL::ORIGIN + normal / distance_to_origin;
  }

  // add facets
  void add_facets(builder &B)
  {
    Vertex_iterator v;
    for(v = m_pMesh->vertices_begin();
      v != m_pMesh->vertices_end();
      v++)
    {
      B.begin_facet();
      HV_circulator he = v->vertex_begin();
      HV_circulator end = he;
      CGAL_For_all(he,end)
	B.add_vertex_to_facet(m_face_map[he->facet()]);
      B.end_facet();
    }
  }

};

template <class Polyhedron, class Kernel>
class Dualizer
{
public:
  typedef typename Polyhedron::HalfedgeDS HalfedgeDS;
  Dualizer() {}
  ~Dualizer() {}

public:
  void run(Polyhedron &input, Polyhedron &output)
  {
    CModifierDual<HalfedgeDS,Polyhedron,Kernel> dualizer(&input);
    output.delegate(dualizer);
  }
};

#endif // _DUALIZER_
