#include <CGAL/Unique_hash_map.h>
#include <fstream>
#include <vector>
#include <list>


struct VERTEX;
struct EDGE;
struct FACET;

typedef std::list<FACET>::iterator FIT;
typedef std::vector<VERTEX>::iterator VIT;
typedef CGAL::Unique_hash_map<VIT, CGAL::Unique_hash_map<VIT, EDGE> > MAP;

struct VERTEX {
  double point[3];
  int index;
  int incident_facets;

  VERTEX(double p[3]) : index(0), incident_facets(0) {
    point[0] = p[0];
    point[1] = p[1];
    point[2] = p[2];
  }
};

struct EDGE {
  std::list<FIT> facets;
};

struct FACET {
  bool erased;
  int points[3];
  int original[3];

  FACET(int p[3]) : erased(false) {
    original[0] = points[0] = p[0];
    original[1] = points[1] = p[1];
    original[2] = points[2] = p[2];

    if(points[1] < points[0]) 
      exchange(points[0], points[1]);
    if(points[2] < points[1]) 
      exchange(points[2], points[1]);
    if(points[1] < points[0]) exchange(points[0], points[1]);    
  }

  void exchange(int& a, int& b) {
    int x = a;
    a = b;
    b = x;
  }
};

int main(int argc, char* argv[]) {

  if(argc != 2)
    return 1;

  std::ifstream trunk(argv[1]);
  
  char buffer[80];
  int nv,nf,nc,nfl;
  std::vector<VERTEX> vertices;
  std::list<FACET> facets;
  //  MAP edges;

  trunk.getline(buffer,80);
  int version = std::atoi(buffer+10);
  trunk.getline(buffer,80);
  trunk >> nv;
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nf;
  if(version > 3) {
    trunk.getline(buffer,80);
    trunk.getline(buffer,80);
    trunk >> nfl;
  }
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);
  trunk >> nc;

  trunk.getline(buffer,80);
  for(int i=0;i<nc;++i)
    trunk.getline(buffer,80);

  trunk.getline(buffer,80);

  std::cerr << "number of vertices " << nv << std::endl;
  for(int i=0;i<nv;++i) {
    double p[3];
    trunk >> p[0] >> p[1] >> p[2];
    vertices.push_back(p);
  }
  
  std::cerr << "number of facets " << nf << std::endl;
  trunk.getline(buffer,80);
  trunk.getline(buffer,80);

  int face[3], col, man;
  for(int i=0;i<nf;++i){
    if(i%1000 == 0)
      std::cerr << i << std::endl;
    trunk >> face[0] >> face[1] >> face[2] >> col >> man;

    VIT vi0 = vertices.begin()+face[0];
    VIT vi1 = vertices.begin()+face[1];
    VIT vi2 = vertices.begin()+face[2];

    ++(vi0->incident_facets);
    ++(vi1->incident_facets);
    ++(vi2->incident_facets);

    FACET f(face);
    facets.push_back(face);
    FIT fi = --facets.end();
    //    edges[vi0][vi1].facets.push_back(fi);
    //    edges[vi0][vi2].facets.push_back(fi);
    //    edges[vi1][vi2].facets.push_back(fi);
  }  

  /*  
  std::list<FIT> border;
  for(FIT fi = facets.begin(); fi != facets.end(); ++fi) {
    VIT vi0 = vertices.begin()+fi->points[0];
    VIT vi1 = vertices.begin()+fi->points[1];
    VIT vi2 = vertices.begin()+fi->points[2];
    if(edges[vi0][vi1].facets.size() < 2 ||
       edges[vi0][vi2].facets.size() < 2 ||
       edges[vi1][vi2].facets.size() < 2)
      border.push_back(fi);
  }
  
  std::cerr << "borders " << border.size() << std::endl;
  */

  /*
  while(border.size() > 0) {
    FIT fi = border.front();
    border.pop_front();
    //    fi->erased = true;
    //    --nf;
    //    --(vi0->incident_facets);
    //    --(vi1->incident_facets);
    //    --(vi2->incident_facets);
  }
  */

  std::cout << "OFF" << std::endl;
  std::cout << vertices.size() << " "
	    << nf << " 0" << std::endl;
  
  std::vector<VERTEX>::iterator vi;
  for(vi = vertices.begin(); vi != vertices.end(); ++vi)
    std::cout << vi->point[0] << " " 
	      << vi->point[1] << " "
	      << vi->point[2] << std::endl;

  std::list<FACET>::iterator fi;
  for(fi = facets.begin(); fi != facets.end(); ++fi)
    if(!fi->erased)
      std::cout << "3 " << fi->original[0] 
		<< " "  << fi->original[1]
		<< " "  << fi->original[2] << std::endl;
  
  return 0;
}
