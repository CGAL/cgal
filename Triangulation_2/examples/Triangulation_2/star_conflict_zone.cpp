#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/spatial_sort.h>
#include <array>
#include <vector>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Delaunay_triangulation_2<K>  Dt2;
typedef Dt2::Edge Edge;
typedef Dt2::Face_handle Face_handle;
typedef Dt2::Face_circulator Face_circulator;
typedef Dt2::Vertex_handle Vertex_handle;

int main( )
{
  Dt2 dt2;

  dt2.insert(Point_2(0,0));
  dt2.insert(Point_2(10,0));
  dt2.insert(Point_2(0,10));

  std::array<Point_2,3> points = { Point_2(2,2), Point_2(1,0), Point_2(2,2) };

  CGAL::spatial_sort(points.begin(), points.end());

  Face_handle hint;

  std::vector<Face_handle> faces;
  std::vector<Edge> edges;

  assert(dt2.dimension() == 2); // precondition of get_conflicts_and_boundary
  for(const Point_2 p : points){
    faces.clear(); // faster than variables in the scope
    edges.clear();
    dt2.get_conflicts_and_boundary(p,
                                   std::back_inserter(faces),
                                   std::back_inserter(edges),
                                   hint);

    if(faces.empty()){
      std::cout << "point " << p << " already in the triangulation" << std::endl;
    }else{
      // Do something with the faces before the insertion

      Vertex_handle vh = dt2.star_hole(p,
                                       edges.begin(), edges.end(),
                                       faces.begin(), faces.end());
      hint = vh->face(); // we could also take any element of faces

      // Do something with the faces after the insertion
      Face_circulator fc = dt2.incident_faces(vh), done(fc);
      do {
        fc++;
      } while (fc != done);
    }
    draw(dt2);
  }
  return 0;
}
