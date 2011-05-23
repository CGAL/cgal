#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <iostream>
#include <list>
#include <fstream>

typedef CGAL::Cartesian<double>                                      Kernel;
typedef Kernel::Point_3                                              Point;
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor		vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator		vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor		edge_descriptor;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator		out_edge_iterator;


using namespace std;

void k_ring(Polyhedron &P, vertex_descriptor v, size_t k, vector<vertex_descriptor> &neigh_vtx)
{
	neigh_vtx.clear();
	neigh_vtx.push_back(v);
	int idx_lv = 0;								// pointing the neighboring vertices on current level
	int idx_lv_end;

	for (size_t lv = 0; lv < k; lv++)
	{
		idx_lv_end = neigh_vtx.size(); 
		for ( ;idx_lv < idx_lv_end; idx_lv++ )
		{
			vertex_descriptor vd = neigh_vtx[idx_lv];
			out_edge_iterator e, e_end;
			for (boost::tie(e,e_end) = boost::out_edges(vd, P); e != e_end; e++)
			{
				vertex_descriptor vt = boost::target(*e, P);
				vector<vertex_descriptor> ::iterator result = find(neigh_vtx.begin(), neigh_vtx.end(), vt);
				if (result == neigh_vtx.end())
				{
					neigh_vtx.push_back(vt);
				}
			}
		}
	}

}


int main() {

	/////////// input OFF file /////////////////////////
	Polyhedron P;
	string filename;
	cout << "Input filename: ";
	cin >> filename;
	ifstream file;
	file.open(&filename[0]);
	file >> P;
	cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << endl;

	////////// extract k-ring vertices ////////////////////////
	vector<vertex_descriptor> neigh_vtx;
	cout << "Determine level of k-ring: ";
	size_t k;
	cin >> k;
	//cout << "Determine vertex index: ";
	k_ring(P, P.vertices_begin(), k, neigh_vtx);

	///////// output indices of neighboting vertices ///////////////////////
	vertex_iterator vb, ve;		 // associate indices to the vertices using the "id()" field of the vertex.
	int index = 0;
	// boost::tie assigns the first and second element of the std::pair
	// returned by boost::vertices to the variables vit and ve
	for(boost::tie(vb,ve)=boost::vertices(P); vb!=ve; ++vb ){
		vertex_descriptor  vd = *vb;
		vd->id() = index++;
	}
	cout << endl << neigh_vtx.size() << " neighboring vertices:" << endl;
	for (int i = 0; i < neigh_vtx.size(); i++)
	{
	cout << neigh_vtx[i]->id() << endl;
	}

	return 0;
}