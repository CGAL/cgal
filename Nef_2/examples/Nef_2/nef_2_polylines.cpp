#include <CGAL/Exact_rational.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>


typedef CGAL::Lazy_exact_nt<CGAL::Exact_rational> FT;
typedef CGAL::Simple_cartesian<FT> Kernel;
typedef CGAL::Bounded_kernel<Kernel> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Point Point;

typedef Nef_polyhedron::Explorer Explorer;
typedef Explorer::Face_const_iterator Face_const_iterator;
typedef Explorer::Hole_const_iterator Hole_const_iterator;
typedef Explorer::Halfedge_around_face_const_circulator Halfedge_around_face_const_circulator;
typedef Explorer::Vertex_const_handle Vertex_const_handle;

int main()
{
  Point r1[3] = { Point(20,15), Point(25,5), Point(30,15) };
  //Point s1[3] = { Point(40,15), Point(40,5) };
  Point t1[3] = { Point(20,10), Point(30,10), Point(25, 15) };

  std::list<std::pair<Point*,Point*> > polylines;
  polylines.push_back(std::make_pair(r1+0, r1+3));
  //polylines.push_back(std::make_pair(s1+0, s1+2));
  polylines.push_back(std::make_pair(t1+0, t1+3));


   Nef_polyhedron RST(polylines.begin(), polylines.end(), Nef_polyhedron::POLYLINES);

   std::cout << RST << std::endl;

    Explorer explorer = RST.explorer();


    // The first face is the infinite one. It has no outer face cycle but only holes

    Face_const_iterator fit = explorer.faces_begin();
    std::cout << "explorer.mark(explorer.faces_begin()) "  << ((explorer.mark(fit))? "is part of polygon" :  "is not part of polygon") << std::endl;
    for(Hole_const_iterator hit = explorer.holes_begin(fit); hit != explorer.holes_end(fit); hit++){
        std::cout << " A hole" << std::endl;
        Halfedge_around_face_const_circulator hafc(hit), done(hit);
        do{
          Vertex_const_handle vh = explorer.target(hafc);
          std::cout << explorer.point(vh) << " [" << explorer.point(vh).x().exact() << " | " << explorer.point(vh).y().exact() << "],  " ;
          hafc++;
        }while(hafc != done);
        std::cout << std::endl;
      }


    // The other faces have outer face cycles, and they may have holes
    for( fit++;
        fit != explorer.faces_end();
        fit++){

      Halfedge_around_face_const_circulator hafc = explorer.face_cycle(fit), done(hafc);
      std::cout << "face: " << ((explorer.mark(fit))? "is part of polygon" :  "is not part of polygon") << std::endl;
      do{
        Vertex_const_handle vh = explorer.target(hafc);
          std::cout << explorer.point(vh) << " [" << explorer.point(vh).x().exact() << " | " << explorer.point(vh).y().exact() << "],  " ;
          hafc++;
      }while(hafc != done);
      std::cout << std::endl;

      for(Hole_const_iterator hit = explorer.holes_begin(fit); hit != explorer.holes_end(fit); hit++){
        std::cout << " A hole" << std::endl;
        Halfedge_around_face_const_circulator hafc(hit), done(hit);
        do{
          Vertex_const_handle vh = explorer.target(hafc);
          std::cout << explorer.point(vh) << " [" << explorer.point(vh).x().exact() << " | " << explorer.point(vh).y().exact() << "],  " ;
          hafc++;
        }while(hafc != done);
        std::cout << std::endl;
      }

    }


  return 0;
}
