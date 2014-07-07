// used by test_intersections_3.cpp

#include <CGAL/Modifier_base.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

const int cube[][3] = { { 0, 1, 3 },
                        { 3, 1, 2 },
                        { 0, 4, 1 },
                        { 1, 4, 5 },
                        { 3, 2, 7 },
                        { 7, 2, 6 },
                        { 4, 0, 3 },
                        { 7, 4, 3 },
                        { 6, 4, 7 },
                        { 6, 5, 4 },
                        { 1, 5, 6 },
                        { 2, 1, 6 } };

template <typename Polyhedron>
struct Build_bbox_mesh :
  public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
  CGAL::Bbox_3 bbox;

public:
  Build_bbox_mesh(CGAL::Bbox_3 b)
    : bbox(b)
  {}

  typedef typename Polyhedron::HalfedgeDS HDS;

  void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( 8, 12, 24);
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;
    B.add_vertex( Point( bbox.xmin(), bbox.ymin(), bbox.zmin())); // -1 -1 -1
    B.add_vertex( Point( bbox.xmin(), bbox.ymax(), bbox.zmin())); // -1 1 -1
    B.add_vertex( Point( bbox.xmax(), bbox.ymax(), bbox.zmin())); // 1 1 -1
    B.add_vertex( Point( bbox.xmax(), bbox.ymin(), bbox.zmin())); // 1 -1 -1
    B.add_vertex( Point( bbox.xmin(), bbox.ymin(), bbox.zmax())); // -1 -1 1
    B.add_vertex( Point( bbox.xmin(), bbox.ymax(), bbox.zmax())); // -1 1 1
    B.add_vertex( Point( bbox.xmax(), bbox.ymax(), bbox.zmax())); // 1 1 1
    B.add_vertex( Point( bbox.xmax(), bbox.ymin(), bbox.zmax())); // 1 -1 1
    for(int i = 0; i < 12; ++i) {
      B.begin_facet();
      B.add_vertex_to_facet( cube[i][0]);
      B.add_vertex_to_facet( cube[i][1]);
      B.add_vertex_to_facet( cube[i][2]);
      B.end_facet();
    }
    B.end_surface();
  }
};

template <typename Polyhedron>
const Polyhedron create_bbox_mesh(const CGAL::Bbox_3& bbox)
{
  Polyhedron result;
  Build_bbox_mesh<Polyhedron> build_bbox(bbox);
  result.delegate(build_bbox);
  return result;
}
