/*! \file
 * A construction of a Topological_map with 2 dimensional points.
 */

#include <CGAL/Cartesian.h>

#include <iostream>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Topological_map.h>

typedef CGAL::Point_2<CGAL::Cartesian<double> >   Point;
typedef CGAL::Segment_2<CGAL::Cartesian<double> > Curve;

typedef CGAL::Pm_dcel<CGAL::Pm_vertex_base<Point>,
                      CGAL::Pm_halfedge_base<Curve>,
                      CGAL::Pm_face_base > Dcel;  

typedef CGAL::Topological_map<Dcel> Tpm;

typedef Tpm::Halfedge_handle Halfedge_handle;
typedef Tpm::Vertex_handle   Vertex_handle;
typedef Tpm::Face_handle     Face_handle;

int main()
{
  Point p1(0, 0), p2(0, 1), p3(1, 1);
  Curve c1(p1, p2), c2(p2, p3), c3(p3, p1);

  Tpm t;

  Face_handle uf = t.unbounded_face();

  Halfedge_handle e1 = t.insert_in_face_interior(uf);
  e1->set_curve(c1);e1->twin()->set_curve(c1);
  e1->source()->set_point(p1);
  e1->target()->set_point(p2);
  CGAL_assertion(t.is_valid());
  std::cout << "Curve " << c1 << " inserted in unbounded face interior"
            << std::endl;

  Halfedge_handle e2 = t.insert_from_vertex(e1);
  e2->set_curve(c2);
  e2->twin()->set_curve(c2);
  e2->target()->set_point(p3);
  CGAL_assertion(t.is_valid());
  std::cout << "Curve " << c2 << " inserted from vertex "
            << e1->target()->point() << std::endl;

  Halfedge_handle e3 = t.insert_at_vertices(e2, e1->twin());
  e3->set_curve(c3);
  e3->twin()->set_curve(c3);
  CGAL_assertion(t.is_valid());
  std::cout << "Curve " << c3 << " inserted between vertices "
            << e2->target()->point() << " and "
            << e1->twin()->target()->point() << std::endl;

  Face_handle nf = e3->face();

  // Print the vertices of the new face. Their order should be 1 1, 0 0, 0 1
  Tpm::Ccb_halfedge_circulator cc = nf->outer_ccb();
  std::cout << std::endl << "Order of vertices in new face : "
            << cc->source()->point();
  for (++cc; cc != nf->outer_ccb(); ++cc)
      std::cout << ", " << cc->source()->point();
  std::cout << std::endl;

  return 0;
}
