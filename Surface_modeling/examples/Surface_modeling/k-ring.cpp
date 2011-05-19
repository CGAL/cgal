#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <fstream>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

typedef Polyhedron::Vertex_const_handle                      Vertex_handle;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;



void load_OFF(Polyhedron &P)
{

}


void extract_k_ring(const Polyhedron &P, Vertex_handle v, int k)
{
  // The map stores the distance of the vertex, 
  // and serves at the same time to find out if we visited the vertex already
  std::map<Vertex_handle, int> dist;
  
  // The queue stores the vertices in increasing distance
  std::list<Vertex_handle> queue;
  dist[v] = 0;
  queue.push_back(v);

  while(! queue.empty()){
    v = queue.front();
    int d = dist[v];
    // When we encounter the first vertex at distance k
    // we can be sure that the points in queue are exactly those at distance k
    if(d == k){
      std::cerr << queue.size() << " vertices at distance " << k << std::endl;
      break;
    }
    queue.pop_front();
    HV_circulator wc = v->vertex_begin(), done(wc);
    do {
      Vertex_handle wh = wc->opposite()->vertex();
      if(dist.find(wh) == dist.end()){
        dist[wh] = d+1;
        queue.push_back(wh);
      }
      ++wc;
    }while(wc != done);
    
  }
}

int main() {
  Polyhedron P;
  std::cin >> P;

  std:: cout << P.size_of_vertices() << " vertices;  " << P.size_of_facets() << "facets" << std::endl;
                                                \
  int k=3;

  extract_k_ring(P, P.vertices_begin(), k);
  return 0;
}

