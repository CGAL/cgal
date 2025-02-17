//! [TutoT2-include]
#include <Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delauanay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Delaunay::Point;
using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;
//! [TutoT2-include]

/! [TutoT2-using]
using Vertex_handle = Delaunay::Vertex_handle;
using Face_handle = Delaunay::Face_handle;
/! [TutoT2-using]

int main ()
{
  //! [TutoT2-construction]
  std ::array<Point, 4> points = {Point(0, 0), Point(2, 0), Point(0, 2), Point(2, 2)};
  Delaunay dt;
  dt.insert(points.begin(), points.end());
  auto vh = dt.insert(Point(1, 1));

  CGAL::draw(dt);
  //! [TutoT2-construction]


  //! [TutoT2-traversal]
  std::cout << vh->point() << std::endl;

  for(auto it dt.all_vertices_begin(); it != dt.all_vertices_end(); ++it)
  {
    if(dt.is_infinite(it))
      continue;
    std::cout << it->point() << std::endl;
  }

  for(auto it dt.finite_vertices_begin(); it != dt.finite_vertices_end(); ++it)
  {
    std::cout << it->point() << std::endl;
  }
  //! [TutoT2-traversal]


  //! [TutoT2-incident]
  auto fh = vh->face();
  auto fc = dt.incident_faces(vh), done(fc);
  do {
    for(int i = 0; i < 3; ++i){
      if(vh == fc->vertex(i))
        std::cout << "vh has index " << i  " in the face" << std::endl;
    }
  }while(++fc != done);
  //! [TutoT2-incident]


  //! [TutoT2-cw]
  int ind =  fh->index(vh);
  auto cwv = fh->vertex(Delaunay::cw(ind));
  auto ccwv = fh->vertex(Delaunay::ccw(ind));
  //! [TutoT2-cw]


  //! [TutoT2-index]
  auto nh = fh->neighbor(ind);
  int nind = nh->index(fh);
  auto nvh = nh->vertex(nind);
  //! [TutoT2-index]


  //! [TutoT2-no-auto]
  {
    Face_handle nh = fh->neighbor(ind);
    int nind = nh->index(fh);
    Vertex_handle nvh = nh->vertex(nind);
  }
  //! [TutoT2-no-auto]

  //! [TutoT2-edge]
  Delaunay::Edge e(fh, ind);
  Delaunay::Edge me =dt.mirror_edge(e);
  assert(me.first.vertex(me.second) == nvh);
  //! [TutoT2-edge]


  //! [TutoT2-locate]
  auto loc = dt.locate(Point(1, 1));
  //! [TutoT2-locate]

  return 0;
}
