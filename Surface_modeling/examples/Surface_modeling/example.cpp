

#include <CGAL/Deform_mesh.h>




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

  CGAL::Deform_mesh deform(P);          
	
	// takes arbitrary vertex as handle
	deform.handles(P.vertices_begin(), P.vertices_begin());

	// determine the k-ring
	deform.region_of_interest(P.vertices_begin(), P.vertices_begin(), 2);

	// does the precomputation
	deform.preprocess();

	// displaces the handle by the Vector_3   v - origin
	Vector translation = P.vertices_begin()->point() - CGAL::ORIGIN;
	deform(P.vertices_begin(), translation);

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
