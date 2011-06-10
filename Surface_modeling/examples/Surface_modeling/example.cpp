#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Taucs_solver_traits.h>
#include <CGAL/Deform_mesh_BGL.h>

typedef CGAL::Cartesian<double>                                                      Kernel;
typedef Kernel::Vector_3                                                             Vector;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>                  Polyhedron;
typedef CGAL::Deform_mesh_BGL<Polyhedron, CGAL::Taucs_solver_traits<double> >        Deform_mesh;

typedef boost::graph_traits<Polyhedron>::vertex_iterator		                         vertex_iterator;


int main() {
	
	// read an off file
	Polyhedron P;           //source mesh
  std::cout << "Please input filename: ";
  std::string filename;
  std::cin >> filename;
	std::string fullname = filename + ".off";
	std::ifstream input(&fullname[0]);
	input >> P;
  std::cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << std::endl;
	input.close();

  Deform_mesh deform(P);          
	
  vertex_iterator vb, ve;
  boost::tie(vb,ve) = boost::vertices(P);
	// takes arbitrary vertex as handle
	deform.handles(vb, vb);

	// determine the k-ring
	deform.region_of_interest(vb, vb, 1);

	// does the precomputation
	deform.preprocess();

	// displaces the handle by the Vector_3   v - origin
	Vector translation = (*vb)->point() - CGAL::ORIGIN;
	deform(vb, translation);

	// write the polyhedron to a file
	std::string outname = filename + "_arap.off";
  std::cout << "output to " << outname << "...";
	std::ofstream output(&outname[0]);
	output << deform.polyhedron;
	output.close();
  std::cout << "done." << std::endl;
  output.close();
	return 0;
}
