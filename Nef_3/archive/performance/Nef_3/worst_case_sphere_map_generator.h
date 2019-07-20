#include <iostream>

template <typename Kernel_>
class worst_case_sphere_map_generator {

  typedef Kernel_ Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Plane_3 Plane_3;

 private:
  std::ostream& out;
  
 public:

  worst_case_sphere_map_generator(std::ostream &o) : out(o) {}

  void print(int n) {

    int x = 2*(n+1);

    out << "Selective Nef Complex" << std::endl;
    out << "standard" << std::endl;
    out << "vertices " << 2*n+3 << std::endl;
    out << "halfedges " << 6*n+6 << std::endl;
    out << "facets " << 2*n+2 << std::endl;
    out << "volumes " << 1 << std::endl;
    out << "shalfedges " << 6*n+6 << std::endl;
    out << "shalfloops " << 0 << std::endl;
    out << "sfaces " << 2*n+3 << std::endl;

    // print vertices

    std::vector<Point_3> points;
    points.push_back(Point_3(0,0,0));
    for(int i=0; i<n+1; ++i) {
      points.push_back(Point_3(n,-n+i*2,-n));
      points.push_back(Point_3(n,-n+i*2, n));
    }

    out << "0 { 0 " << 2*n+1 << ", 0 " << 2*n+1 << ", 0 0, -2 | "
	<< *points.begin() << " } 1" << std::endl;
    for(int i=0; i<2*(n+1); ++i) {
      out << i+1 << " { " 
	  << 2*(n+i)+2 << " " << 2*(n+i)+3 << ", "
	  << 2*(n+i)+2 << " " << 2*(n+i)+3 << ", "
	  << i+1 << " " << i+1 << ", -2 | "
	  << points[i+1] << " } 1" << std::endl;
    }
    
    // print halfedges

    for(int i=0; i<n+1;++i) {
      out << 2*i   << " { " << x+4*i << ", " 
	  << 0 << ", 0 "
	  << 2*i << " | " 
	  << points[1+2*i] - points[0] << " } 1" << std::endl;

      out << 2*i+1 << " { " << x+4*i+2 << ", " 
	  << 0 << ", 0 "
	  << 2*i+1 << " | " 
	  << points[2+2*i] - points[0] << " } 1" << std::endl;
    }

    for(int i=0; i<n+1;++i) {
      out << x+4*i << " { " << 2*i << ", " 
	  << 1+2*i << ", 0 "
	  << x+4*i+1 << " | " 
	  << points[0] - points[1+2*i] << " } 1" << std::endl;

      out << x+4*i+1 << " { " << x+4*i+3 << ", " 
	  << 1+2*i << ", 0 "
	  << x+4*i << " | " 
	  << points[2+2*i] - points[1+2*i] << " } 1" << std::endl;

      out << x+4*i+2 << " { " << 2*i+1 << ", " 
	  << 2+2*i << ", 0 "
	  << x+4*i+2 << " | " 
	  << points[0] - points[2+2*i] << " } 1" << std::endl;

      out << x+4*i+3 << " { " << x+4*i+1 << ", " 
	  << 2+2*i << ", 0 "
	  << x+4*i+3 << " | " 
	  << points[1+2*i] - points[2+2*i] << " } 1" << std::endl;



    }

    // print facets

    for(int i=0; i<n+1; ++i) {
      out << 2*i << " { " << 2*i+1 << ", "
	  << 2*i << " , , 0 | "
	  << Plane_3(points[2*i+1], points[0], points[2*i+2]) << " } 1" << std::endl;

      out << 2*i+1 << " { " << 2*i << ", "
	  << 2*i+1 << " , , 0 | "
	  << Plane_3(points[2*i+2], points[0], points[2*i+1]) << " } 1" << std::endl;
    }

    // print volume

    out << "0 { 0 } 0" << std::endl;

    // print sedges

    for(int i=0; i<n+1; ++i) {
      out << 2*i << " { " << 2*i+1 << ", "
	  << 2*i+1 << ", " << 2*i+1 << ", "
	  << 2*i   << ", " << 0 << ", "
	  << x+4*i << ", " << x+4*i+2 << ", " 
	  << 2*i << " | " 
	  << Plane_3(points[2*i+2], points[0], points[2*i+1]) << " } 1" << std::endl;

      out << 2*i+1 << " { " << 2*i << ", "
	  << 2*i << ", " << 2*i << ", "
	  << 2*i+1 << ", " << 0 << ", "
	  << x+4*i+3 << ", " << x+4*i+1 << ", " 
	  << 2*i+1 << " | " 
	  << Plane_3(points[2*i+1], points[0], points[2*i+2]) << " } 1" << std::endl;
    }      

    for(int i=0; i<n+1; ++i) {
      out << x+4*i << " { " << x+4*i+1 << ", "
	  << x+4*i+1 << ", " << x+4*i+1 << ", "
	  << x+4*i+1 << ", " << 2*i+1 << ", "
	  << x+4*i+2 << ", " << 2*i << ", " 
	  << 2*i << " | " 
	  << Plane_3(points[2*i+2], points[0], points[2*i+1]) << " } 1" << std::endl;

      out << x+4*i+1 << " { " << x+4*i << ", "
	  << x+4*i << ", " << x+4*i << ", "
	  << x+4*i << ", " << 2*i+1 << ", "
	  << 2*i+1 << ", " << x+4*i+3 << ", " 
	  << 2*i+1 << " | " 
	  << Plane_3(points[2*i+1], points[0], points[2*i+2]) << " } 1" << std::endl;

      out << x+4*i+2 << " { " << x+4*i+3 << ", "
	  << x+4*i+3 << ", " << x+4*i+3 << ", "
	  << x+4*i+2 << ", " << 2*i+2 << ", "
	  << 2*i << ", " << x+4*i << ", " 
	  << 2*i << " | " 
	  << Plane_3(points[2*i+2], points[0], points[2*i+1]) << " } 1" << std::endl;

      out << x+4*i+3 << " { " << x+4*i+2 << ", "
	  << x+4*i+2 << ", " << x+4*i+2 << ", "
	  << x+4*i+3 << ", " << 2*i+2 << ", "
	  << x+4*i+1 << ", " << 2*i+1 << ", " 
	  << 2*i+1 << " | " 
	  << Plane_3(points[2*i+1], points[0], points[2*i+2]) << " } 1" << std::endl;
    }

    // print sfaces

    out << "0 { 0, ";
    for(int i=0; i<n+1; ++i)
      out << 2*i << " ";
    out << ", , , 0 } 0" << std::endl;
    
    for(int i=0; i<n+1; ++i) {
      out << 2*i+1 << " { " << 2*i+1 << ", "
	  << x+4*i << " , , , 0 } 0" << std::endl;
      
      out << 2*i+2 << " { " << 2*i+2 << ", "
	  << x+4*i+2 << " , , , 0 } 0" << std::endl;
    }
    
    out << "/* end Selective Nef complex */" << std::endl;
  }
};
