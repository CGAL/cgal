

#include <CGAL/Deform_mesh.h>




int main() {
	
	// read an off file
	Polyhedron P;           //source mesh
	cout << "Please input filename: ";
	string filename;
	cin >> filename;
	string fullname = filename + ".off";
	ifstream input(&fullname[0]);
	input >> P;
	cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << endl;
	input.close();

	Deform_mesh deform(P);          
	
	// takes arbitrary vertex as handle
	deform.handles(P.vertices_begin(), P.vertices_begin());

	// determine the k-ring
	deform.region_of_interest(P.vertices_begin(), P.vertices_begin(), 2);

	// does the precomputation
	deform.preprocess();

	// displaces the handle by the Vector_3   v - origin
	Vector translation = P.vertices_begin()->point() - CGAL::ORIGIN;
	deform(P.vertices_begin()->point(), translation);

	// write the polyhedron to a file
	string outname = filename + "_arap.off";
	cout << "output to " << outname << "...";
	ofstream output(&outname[0]);
	output << P;
	output.close();
	cout << "done." << endl;
	return 0;
}
