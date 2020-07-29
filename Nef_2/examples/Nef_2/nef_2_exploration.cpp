#include <CGAL/Exact_integer.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>

typedef CGAL::Exact_integer RT;
typedef CGAL::Filtered_extended_homogeneous<RT> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;

typedef Nef_polyhedron::Point Point;  // a standard point
typedef Nef_polyhedron::Line  Line;
typedef Nef_polyhedron::Explorer Explorer;

typedef Explorer::Face_const_iterator Face_const_iterator;
typedef Explorer::Hole_const_iterator Hole_const_iterator;
typedef Explorer::Halfedge_around_face_const_circulator Halfedge_around_face_const_circulator;
typedef Explorer::Vertex_const_handle Vertex_const_handle;


void explore(std::string s, const Nef_polyhedron&  poly)
{
  std::cout << "Explore: " << s << std::endl;

  Explorer explorer = poly.explorer();
  int i = 0;
  for(Face_const_iterator fit = explorer.faces_begin(); fit != explorer.faces_end(); ++fit, i++){
    std::cout << "\nFace " << i << std::endl
              << ( explorer.mark(fit)? "* is" : "* is not") << " marked" << std::endl;

    // explore the outer face cycle if it exists
    Halfedge_around_face_const_circulator hafc = explorer.face_cycle(fit);
    if(hafc == Halfedge_around_face_const_circulator()){
      std::cout << "* has no outer face cycle" << std::endl;
    } else {
      std::cout << "* outer face cycle:\n";
      std::cout << "  - halfedges around the face: ";
      Halfedge_around_face_const_circulator done(hafc);
      do {
        char c = (explorer.is_frame_edge(hafc))?'f':'e';
        std::cout << c;
        ++hafc;
      }while (hafc != done);
      std::cout << " ( f = frame edge, e = ordinary edge)" << std::endl;

      std::cout << "  - vertices around the face:\n";
      do {
        Vertex_const_handle vh = explorer.target(hafc);
        if (explorer.is_standard(vh)){
          std::cout << "      " << explorer.point(vh) << std::endl;
        }else{
          std::cout << "      " << explorer.ray(vh) << std::endl;
        }
        ++hafc;
      }while (hafc != done);
    }

    // explore the holes if the face has holes
    Hole_const_iterator hit = explorer.holes_begin(fit), end = explorer.holes_end(fit);
    if(hit == end){
      std::cout << "* has no hole" << std::endl;
    }else{
      std::cout << "* has holes" << std::endl;
      for(; hit != end; hit++){
        Halfedge_around_face_const_circulator hafc(hit), done(hit);
        std::cout << "  - halfedges around the hole: ";
        do {
          char c = (explorer.is_frame_edge(hafc))?'f':'e';
          std::cout << c;
          ++hafc;
        }while (hafc != done);
        std::cout << " ( f = frame edge, e = ordinary edge)" << std::endl;
      }
    }
  }
  std::cout << "done\n" << std::endl;
}

int main() {

  CGAL::set_pretty_mode(std::cout);
  Nef_polyhedron N0(Nef_polyhedron::COMPLETE);
  explore("complete", N0);

  Line l(0,1,-2); // l : 0x + y - 2 = 0
  Nef_polyhedron N1(l,Nef_polyhedron::INCLUDED);
  explore("line", N1);

  Point p1(1,1), p2(10,1), p3(10,10);
  Point triangle[3] = { p1, p2, p3 };
  Nef_polyhedron N2(triangle, triangle+3);
  explore("triangle", N2);

  {
    Point p1(4,2), p2(6,2), p3(6,4);
    Point triangle[3] = { p1, p2, p3 };
    Nef_polyhedron N3(triangle, triangle+3);
    N2 -= N3;
  }
  {
    Point p1(7,2), p2(9,2), p3(9,6);
    Point triangle[3] = { p1, p2, p3 };
    Nef_polyhedron N3(triangle, triangle+3);
    N2 -= N3;
  }

  explore("triangle with two triangular holes", N2);

  return 0;
}
