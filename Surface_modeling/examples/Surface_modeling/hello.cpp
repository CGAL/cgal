#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <string>
#include <fstream>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Polyhedron::HalfedgeDS              HalfedgeDS;

typedef Polyhedron::Vertex_const_handle                      Vertex_handle;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;

using namespace std;

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
	Build_triangle() {}
	void operator()( HDS& hds) {

		string filename;
		cin >> filename;
		ifstream file;
		file.open(&filename[0]);
		string title;
		file >> title;
		if (title != "OFF")
		{
			cout << "Wrong format!" << endl;
			return;
		}
		int nV, nF, nE;
		file >> nV >> nF >> nE;
		// Postcondition: `hds' is a valid polyhedral surface.
		CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
		B.begin_surface( nV, nF);
		typedef typename HDS::Vertex   Vertex;
		typedef typename Vertex::Point Point;
		for (int idx = 0; idx < nV; idx ++)
		{
			double x, y, z;
			file >> x >> y >> z;
			B.add_vertex( Point( x, y, z));
		}
		for (int fidx = 0; fidx < nF; fidx ++)
		{
			int nP;
			file >> nP;
			B.begin_facet();
			for (int i = 0; i < nP; i ++)
			{
				int idx;
				file >> idx;
				B.add_vertex_to_facet(idx);
			}
			B.end_facet();
		}
		B.end_surface();
		file.close();
		
	}
};


void LoadOFF(Polyhedron &P)
{
	Build_triangle<HalfedgeDS> triangle;
	P.delegate( triangle);
}


void k_ring(Polyhedron &P, Vertex_handle v, size_t k, vector<Vertex_handle> &neigh_vtx)
{
	neigh_vtx.clear();
	neigh_vtx.push_back(v);
	int idx_lv = 0;    // pointing the neighboring vertices on current level
	int idx_lv_end;


	for (size_t lv = 0; lv < k; lv++)
	{
		idx_lv_end = neigh_vtx.size(); 
		for ( ;idx_lv < idx_lv_end; idx_lv++ )
		{
			Vertex_handle vh = neigh_vtx[idx_lv];
			HV_circulator wc = vh->vertex_begin(), done(wc);
			do {
				Vertex_handle wh = wc->opposite()->vertex();
				vector<Vertex_handle> ::iterator result = find(neigh_vtx.begin(), neigh_vtx.end(), wh);
				if (result == neigh_vtx.end())
				{
					neigh_vtx.push_back(wh);
				}
				++wc;
			}while(wc != done);
		}
	}
	
}

int main() {
	Polyhedron P;
	string filename;
	cin >> filename;
	ifstream file;
	file.open(&filename[0]);
	file >> P;
	cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << endl;

	vector<Vertex_handle> neigh_vtx;
	cout << "Determine level of k-ring: ";
	size_t k;
	cin >> k;
	cout << "Determine vertex index: ";
	k_ring(P, P.vertices_begin(), k, neigh_vtx);
	cout << endl << neigh_vtx.size() << " neighboring vertices:" << endl;
	/*for (int i = 0; i < neigh_vtx.size(); i++)
	{
		cout << neigh_vtx[i] << endl;
	}*/

	return 0;
}

