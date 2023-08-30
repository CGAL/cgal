#include <CGAL/Exact_rational.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Polygon_2.h>

typedef CGAL::Lazy_exact_nt<CGAL::Exact_rational> FT;
typedef CGAL::Simple_cartesian<FT> Kernel;
typedef CGAL::Bounded_kernel<Kernel> Bounded_kernel;
typedef CGAL::Nef_polyhedron_2<Bounded_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point Point;

typedef Nef_polyhedron::Explorer Explorer;
typedef Explorer::Face_const_iterator Face_const_iterator;
typedef Explorer::Hole_const_iterator Hole_const_iterator;
typedef Explorer::Halfedge_around_face_const_circulator Halfedge_around_face_const_circulator;
typedef Explorer::Vertex_const_handle Vertex_const_handle;
typedef Explorer::Isolated_vertex_const_iterator Isolated_vertex_const_iterator;

void print(Explorer& explorer, Vertex_const_handle vh, bool exact)
{
    std::cout << explorer.point(vh);
    if (exact)
        std::cout << " [" << explorer.point(vh).x().exact() << " | " << explorer.point(vh).y().exact() << "]";
    std::cout << ", ";
}

void print(const Nef_polyhedron& poly, bool exact = false)
{


    Explorer explorer = poly.explorer();

    // The first face is the infinite one. It has no outer face cycle but only holes
    Face_const_iterator fit = explorer.faces_begin();
    std::cout << "explorer.mark(explorer.faces_begin()) " << ((explorer.mark(fit)) ? "is part of polygon" : "is not part of polygon") << std::endl;
    for (Hole_const_iterator hit = explorer.holes_begin(fit); hit != explorer.holes_end(fit); hit++) {
        std::cout << " A hole" << std::endl;
        Halfedge_around_face_const_circulator hafc(hit), done(hit);
        do {
            Vertex_const_handle vh = explorer.target(hafc);
            print(explorer, vh, exact);
            hafc++;
        } while (hafc != done);
        std::cout << std::endl;
    }

    for (Isolated_vertex_const_iterator it = explorer.isolated_vertices_begin(fit); it != explorer.isolated_vertices_end(fit); ++it) {
        std::cout << "isolated vertex A" << explorer.point(it) << std::endl;
    }

    // The other faces have outer face cycles, and they may have holes
    for (fit++;
        fit != explorer.faces_end();
        fit++) {
        for (Isolated_vertex_const_iterator it = explorer.isolated_vertices_begin(fit); it != explorer.isolated_vertices_end(fit); ++it) {
            std::cout << "isolated vertex B" << std::endl;
        }
        Halfedge_around_face_const_circulator hafc = explorer.face_cycle(fit), done(hafc);
        std::cout << "face: " << ((explorer.mark(fit)) ? "is part of polygon" : "is not part of polygon") << std::endl;
        do {
            Vertex_const_handle vh = explorer.target(hafc);
            print(explorer, vh, exact);
            hafc++;
        } while (hafc != done);
        std::cout << std::endl;
        for (Hole_const_iterator hit = explorer.holes_begin(fit); hit != explorer.holes_end(fit); hit++) {
            std::cout << " A hole" << std::endl;
            Halfedge_around_face_const_circulator hafc(hit), done(hit);
            do {
                Vertex_const_handle vh = explorer.target(hafc);
                print(explorer, vh, exact);
                hafc++;
            } while (hafc != done);
            std::cout << std::endl;
        }
    }

}



int main()
{
  Point p1[] = {
    Point(0,0),
    Point(100,100),
    Point(0,100)
  };

  Point p2[] = {
    Point(100, 100),
    Point(110, 100),
    Point(100, 110)
  };


  std::list<std::pair<Point*, Point*> > polygons1;
  polygons1.push_back(std::make_pair(p1, p1 + sizeof(p1) / sizeof(Point)));

  Nef_polyhedron poly1(polygons1.begin(), polygons1.end(), Nef_polyhedron::POLYGONS);
  poly1.explorer().check_integrity_and_topological_planarity();

  Nef_polyhedron poly2(p2, p2 + sizeof(p2) / sizeof(Point));
  poly2.explorer().check_integrity_and_topological_planarity();

  Nef_polyhedron intersect = poly1.intersection(poly2);
  intersect.explorer().check_integrity_and_topological_planarity();

  print(intersect);

  Nef_polyhedron comp = intersect.complement();  //leads to crash/exception since topological plane map of intersect is not correct!
  comp.explorer().check_integrity_and_topological_planarity();
  return 0;
}
