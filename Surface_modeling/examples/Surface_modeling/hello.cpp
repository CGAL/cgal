#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <string>
#include <fstream>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;
typedef Polyhedron::HalfedgeDS              HalfedgeDS;
typedef Kernel::Vector_3                                     Vector;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

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

void Extract_OneRing(Polyhedron &P, vector<vector<int > > &neigh_vtx)
{
	neigh_vtx.clear();
	neigh_vtx.resize(P.size_of_vertices());
	for (int i = 0; i <neigh_vtx.size(); i++)
	{
		neigh_vtx[i].push_back(i);
	}
	map<Vertex_iterator, int> vidx;
	int idx = 0;
	for (Vertex_iterator vit = P.vertices_begin(); vit != P.vertices_end(); vit++)
	{
		vidx[vit] = idx++;
	}

	vector<map<int, int> > edges(P.size_of_vertices());
	for (Facet_iterator fit = P.facets_begin(); fit != P.facets_end(); fit++)
	{
		HF_circulator h = fit->facet_begin();
		HF_circulator h_next = fit->facet_begin();
		h_next++;
		vector<int> test;
		do  
		{
			if (edges[vidx[h->vertex()]][vidx[h_next->vertex()]] == 0)
			{
				neigh_vtx[vidx[h->vertex()]].push_back(vidx[h_next->vertex()]);
				neigh_vtx[vidx[h_next->vertex()]].push_back(vidx[h->vertex()]);
				edges[vidx[h->vertex()]][vidx[h_next->vertex()]] = 1;
				edges[vidx[h_next->vertex()]][vidx[h->vertex()]] = 1;
			}
			h_next++;
		} while ( ++h != fit->facet_begin() );
	}

}

void k_Ring(Polyhedron &P, size_t k, int idx, vector<int> &neigh_vtx_lv_k)
{


	vector<vector<int > > neigh_vtx;
	Extract_OneRing(P, neigh_vtx);
	neigh_vtx_lv_k = neigh_vtx[idx];          // setting one-ring as the lowest level

	vector<int> border_vtx = neigh_vtx[idx];  // "boundary" vertices on relative levels

	for (size_t level = 1; level < k; level++)
	{
		vector<int> new_border_vtx;
		while (!border_vtx.empty())
		{
			int border_idx = border_vtx.back();
			border_vtx.pop_back();
			for (int i = 0; i < neigh_vtx[border_idx].size(); i++)
			{
				int border_neigh_idx = neigh_vtx[border_idx][i];
				vector<int> ::iterator result = find(neigh_vtx_lv_k.begin(), neigh_vtx_lv_k.end(), border_neigh_idx);
				if (result == neigh_vtx_lv_k.end())   // justify whether the neighboring vertex of border vertex is already included
				{
					neigh_vtx_lv_k.push_back(border_neigh_idx);
					new_border_vtx.push_back(border_neigh_idx);
				}
			}

		}
		border_vtx = new_border_vtx;
	}
}

int main() {
	Polyhedron P;
	LoadOFF(P);
	cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << endl;

	vector<int > neigh_vtx_lv;
	cout << "Determine level of k-ring: ";
	size_t k;
	cin >> k;
	cout << "Determine vertex index: ";
	int idx;
	cin >> idx;
	k_Ring(P, k, idx, neigh_vtx_lv);
	cout << endl << neigh_vtx_lv.size() << " neighboring vertices:" << endl;
	for (int i = 0; i < neigh_vtx_lv.size(); i++)
	{
		cout << neigh_vtx_lv[i] << endl;
	}

	return 0;
}

