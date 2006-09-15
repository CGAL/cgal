#define CGAL_CHECK_EXPENSIVE
#define CGAL_CHECK_EXACTNESS

#include <CGAL/IO/Qt_debug_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_arrangement.h>
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <iostream>
#include <sstream>
#include <fstream>


struct Do_work {
  Do_work(double z, char *name): z_(z), fname_(name){}
  void operator()(Qt_examiner_viewer_2 *q){
    CGAL::Bbox_3 box(std::numeric_limits<double>::max(),
		     std::numeric_limits<double>::max(),
		     std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max(),
		     -std::numeric_limits<double>::max());

    typedef Arrangement_of_spheres_traits_3::Geometric_traits K;
    typedef CGAL::Simple_cartesian<double> DK;
    std::vector<K::Sphere_3> spheres;
    std::ifstream in(fname_);
    while (true){
      char buf[1000];
      in.getline(buf, 1000);
      if (!in) break;
      std::istringstream iss(buf);
      DK::Sphere_3 s;
      iss >> s;
      if (!iss) {
	std::cerr << "Can't parse line " << buf << std::endl;
      } else {
	spheres.push_back(K::Sphere_3(K::Point_3(s.center().x(), 
						 s.center().y(),
						 s.center().z()), 
				      s.squared_radius()));
	box= box+ spheres.back().bbox();
      }
      //std::cout << spheres.back() << std::endl;
    }


    std::cout << "Read " << spheres.size() << " spheres." << std::endl;
    std::cout << "Bounding box is from " << box.zmin() << " to " << box.zmax()
	      << std::endl;
 
    Geometric_traits::FT z= z_;
   
    //
    Slice::T tr(spheres.begin(), spheres.end());
    Slice slice(tr);
    slice.initialize_at(z);
  
    *q << Layer(0);
    slice.draw_rz(q, z);
    q->show_everything();
    //q->redraw();
 

    
    while (true) {
      std::cout << "Enter coordinates: " << std::flush;
      char buf[1000];
      std::cin.getline(buf,1000);
      if (buf[0]== '\0') {
	std::cout << "bye." << std::endl;
	break;
      }
      slice.marked_faces_clear();
      slice.marked_edges_clear();
      slice.marked_vertices_clear();
      double x,y;
      std::istringstream iss(buf);
      iss >> x >> y;
      if (!iss) {
	std::cerr << "Can't parse line." << std::endl;
      } else {
	*q << Layer(1);
	*q << CGAL::RED;
	*q << K::Point_2(x,y);
	//q->redraw();
	Arrangement_of_spheres_traits_3::Sphere_point_3 sp(K::Point_3(x,y,z), 
							  K::Line_3(K::Point_3(x,y,z),
								    K::Vector_3(0,0,1)));
	try {
	  Slice::Face_const_handle f=slice.locate_point(sp);
	  //slice.new_marked_face(f);
	  std::cout << "Found face ";
	  slice.write(f, std::cout) << std::endl;
	  slice.new_marked_face(f);

	} catch (Slice::On_edge_exception e) {
	  std::cout << "On edge!" <<std::endl;
	  slice.new_marked_edge(e.halfedge_handle());
	} catch (Slice::On_vertex_exception v) {
	  std::cout << "On vertex!" <<std::endl;
	  slice.new_marked_vertex(v.vertex_handle());
	}
	//q->clear();
	//slice.draw_rz(q, z);
	*q << Layer(1);
	slice.draw_marked_rz(q,z);
	*q << CGAL::RED;
	*q << K::Point_2(x,y);
	//q->redraw();
      }
    }
  }
  double z_;
  char *fname_;
};


int main(int argc, char *argv[]){
  
  Do_work dw(atof(argv[1]), argv[2]);

  Qt_debug_viewer_2<Do_work> qtd(dw, argc, argv);
  
  return qtd();
}
