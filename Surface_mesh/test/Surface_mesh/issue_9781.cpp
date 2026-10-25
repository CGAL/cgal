   #include <CGAL/Surface_mesh/Surface_mesh.h>
   #include <CGAL/Simple_cartesian.h>
   #include <iostream>

   typedef CGAL::Simple_cartesian<double> K;
   typedef K::Point_3 Point;
   typedef CGAL::Surface_mesh<Point> Mesh;

   int main() {
     Mesh src, dst;
     Mesh::Vertex_index v0 = src.add_vertex(Point(0, 1, 2));
     Mesh::Vertex_index v1 = src.add_vertex(Point(3, 4, 5));

     dst = std::move(src);
     assert(src.is_valid());
     assert(dst.is_valid());
     assert(src.is_empty());
     std::cout << "after move: src.n=" << src.number_of_vertices()<< " src.is_valid=" << src.is_valid() << '\n';

     src.add_vertex(Point(6, 7, 8));

     boost::graph_traits<Mesh>::vertex_descriptor v = *vertices(dst).begin();
     std::cout << dst.point(v) << " (expect 0 1 2)" << std::endl;
     assert(dst.point(v) == Point(0, 1, 2));

     Mesh dst2(std::move(src));
     assert(src.is_valid());
     assert(dst2.is_valid());
     assert(src.is_empty());
     v = *vertices(dst2).begin();
     std::cout << dst2.point(v) << " (expect 6 7 8)" << std::endl;
     return 0;
   }
