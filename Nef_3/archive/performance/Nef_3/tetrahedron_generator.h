#include <CGAL/point_generators_3.h>
#include <list>
#include <vector>

template <typename Kernel_>
class tetrahedron_generator {

  typedef Kernel_ Kernel;
  typedef typename Kernel::RT RT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;
  typedef typename CGAL::Random_points_in_cube_3<Point_3> Point_source;
  
 public:
  static const int sprev[6];
  static const int snext[6];
  static const int source[6];
  static const int sface[12];
  static const int fcycle[8];
  static const int facet[24];
  static const int prev[24];
  static const int next[24];

 private:
  std::ostream& out;
  RT s;
  CGAL::Random r;
  Point_source P;
  std::list<Point_3> points;

  void transform(Point_3& p, int sx, int sy, int sz) {
    p = p.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(sx*2*s+s,sy*2*s+s,sz*2*s+s,2)));
  }
  
  void create_tetra(int sx, int sy, int sz) {
    Point_3 ps[4];

    for(int i=0; i<4; ++i) {
      Point_3 p(*P++);
      transform(p,sx,sy,sz);
      ps[i] = p;
    }
      
    while(ps[0]==ps[1]) {
      ps[1] = *P++;
      transform(ps[1],sx,sy,sz);
    }
      
    if(CGAL::lexicographically_xyz_smaller(ps[1],ps[0])) {
      Point_3 p = ps[0];
      ps[0] = ps[1];
      ps[1] = p;
    }
    
    while(CGAL::collinear(ps[0],ps[1],ps[2])) {
      ps[2] = *P++;
      transform(ps[2],sx,sy,sz);
    }
    
    if(CGAL::lexicographically_xyz_smaller(ps[2],ps[1])) {
      Point_3 p = ps[1];
      ps[1] = ps[2];
      ps[2] = p;
      if(CGAL::lexicographically_xyz_smaller(ps[1],ps[0])) {
	Point_3 p = ps[0];
	ps[0] = ps[1];
	ps[1] = p;
      }       
    }
    
    while(CGAL::orientation(ps[0],ps[1],ps[2],ps[3])==CGAL::COPLANAR) {
      ps[3] = *P++;
      transform(ps[3],sx,sy,sz);
    }
    
    if(CGAL::lexicographically_xyz_smaller(ps[3],ps[1])) {
      Point_3 p = ps[1];
      ps[1] = ps[3];
      ps[3] = p;
      if(CGAL::lexicographically_xyz_smaller(ps[1],ps[0])) {
	Point_3 p = ps[0];
	ps[0] = ps[1];
	ps[1] = p;
      }       
    }
    
    if(CGAL::orientation(ps[0],ps[1],ps[2],ps[3])!=CGAL::POSITIVE) {
      Point_3 p = ps[2];
      ps[2] = ps[3];
      ps[3] = p;
    }
    CGAL_assertion(CGAL::orientation(ps[0],ps[1],ps[2],ps[3])==CGAL::POSITIVE);
    CGAL_assertion(CGAL::orientation(ps[2],ps[0],ps[1],ps[3])==CGAL::POSITIVE);
    CGAL_assertion(CGAL::orientation(ps[3],ps[0],ps[2],ps[1])==CGAL::POSITIVE);
    CGAL_assertion(CGAL::orientation(ps[1],ps[0],ps[3],ps[2])==CGAL::POSITIVE);
    CGAL_assertion(CGAL::orientation(ps[2],ps[1],ps[3],ps[0])==CGAL::POSITIVE);

    for(int i=0; i<4; ++i)
      points.push_back(ps[i]);
  }
  
  void print() {
    CGAL_assertion(points.size()%4 == 0);
    
    std::vector<Vector_3> vectors;
    
    out << "Selective Nef Complex" << std::endl;
    out << "standard" << std::endl;
    out << "vertices " << points.size() << std::endl;
    out << "halfedges " << points.size() * 3 << std::endl;
    out << "facets " << points.size() * 2 << std::endl;
    out << "volumes " << points.size() / 4 + 1 << std::endl;
    out << "shalfedges " << points.size() * 6 << std::endl;
    out << "shalfloops 0" << std::endl;
    out << "sfaces " << points.size() * 2 << std::endl;
    
    // print vertices
    unsigned int idx(0);
    typename std::list<Point_3>::const_iterator pit(points.begin());
    while(pit != points.end()) {
      out << idx << " { " 
		<< idx*3 << " " << idx*3+2 << ", "
		<< idx*6 << " " << idx*6+5 << ", "
		<< idx*2 << " " << idx*2+1 << ", "
		<< "-2 | " 
		<< *pit++ << " } 1" << std::endl;
      ++idx;
    }
    
    // print halfedges
    idx = 0;
    pit = points.begin();
    while(pit != points.end()) {
      Point_3 p[4];
      for(int i=0; i<4; ++i) 
	p[i] = *pit++;
      for(int i=0; i<4; ++i)
	for(int j=0; j<3; ++j) {
	  int opp=(i+j+1)%4;
	  out << idx+i*3+j << " { "
		    << idx+opp*3+((i+3-opp)%4) << ", "
		    << idx/3+i << ", 0 "
		    << idx*2+i*6+(j==2?3:j) << " | " 
		    << p[opp]-p[i] 
		    << " } 1" << std::endl;
	  vectors.push_back(p[opp]-p[i]);
	}
      idx+=12;
    }
    
    // print halffacets
    pit = points.begin();
    for(idx=0;idx<points.size()/4;++idx) {
      out << 8*idx << " { "
		<< 8*idx+1 << ", "
		<< idx*24 << " , , 0 | "
		<< Plane_3(*pit,
			   *pit + vectors[12*idx+1],
			   *pit + vectors[12*idx])
		<< " } 1" << std::endl;
      out << 8*idx+1 << " { "
		<< 8*idx << ", "
		<< idx*24+1 << " , , " << idx+1 << " | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx],
			   *pit + vectors[12*idx+1])
		<< " } 1" << std::endl;    
      
      out << 8*idx+2 << " { "
		<< 8*idx+3 << ", "
		<< idx*24+4 << " , , 0 | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+2],
			   *pit + vectors[12*idx+1])
		<< " } 1" << std::endl;
      out << 8*idx+3 << " { "
		<< 8*idx+2 << ", "
		<< idx*24+5 << " , , " << idx+1 << " | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+1],
			   *pit + vectors[12*idx+2])
		<< " } 1" << std::endl;    
      
      out << 8*idx+4 << " { "
		<< 8*idx+5 << ", "
		<< idx*24+3 << " , , 0 | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+0],
			   *pit + vectors[12*idx+2])
		<< " } 1" << std::endl;
      out << 8*idx+5 << " { "
		<< 8*idx+4 << ", "
		<< idx*24+2 << " , , " << idx+1 << " | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+2],
			   *pit + vectors[12*idx+0])
		<< " } 1" << std::endl;    
      
      ++pit;

      out << 8*idx+6 << " { "
		<< 8*idx+7 << ", "
		<< idx*24+7 << " , , 0 | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+3],
			   *pit + vectors[12*idx+4])
		<< " } 1" << std::endl;
      out << 8*idx+7 << " { "
		<< 8*idx+6 << ", "
		<< idx*24+6 << " , , " << idx+1 << " | "
		<< Plane_3(*pit, 
			   *pit + vectors[12*idx+4],
			   *pit + vectors[12*idx+3])
		<< " } 1" << std::endl;    
      ++pit; ++pit; ++pit;
    }
    
    // print volumes
    out << 0 << " { ";
    for(idx = 0; idx < points.size()/4; ++idx)
      out << idx*8+1 << " ";
    out << "} 0" << std::endl;
    for(idx = 0; idx < points.size()/4; ++idx) 
      out << idx+1 << " { " << idx*8 << " } 1" << std::endl;
    
    // print shalfedges
    idx=0;
    while(idx < points.size()*6) {
      for(int i=0; i<6; ++i) {
	int twin(i%2==0?1:-1);
	out << idx+i << " { " 
		  << idx+i+twin << ", "
		  << idx+sprev[i] << ", "
		  << idx+snext[i] << ", "
		  << idx/2+source[i] << ", "
		  << idx/3+sface[(idx+i)%12] << ", "
		  << idx/24*24+prev[(idx+i)%24] << ", "
		  << idx/24*24+next[(idx+i)%24] << ", "
		  << idx/24*8+facet[(idx+i)%24] << " | "
		  << Plane_3(Point_3(0,0,0), 
			     CGAL::ORIGIN + vectors[idx/2+source[i]],
			     CGAL::ORIGIN + vectors[idx/2+source[i+twin]]) 
		  << " } 1" << std::endl;
      }
      idx+=6;
    }
    
    // print sfaces
    for(idx=0;idx < points.size()/2; ++idx) {
      out << 4*idx << " { " 
		<< idx*2 << ", " 
		<< idx*12 << " , , , "
		<< idx/2+1 << " } 1" << std::endl;
      out << 4*idx+1 << " { " 
		<< idx*2 << ", " 
		<< idx*12+2 << " , , , "
		<< "0 } 0" << std::endl;
      out << 4*idx << " { " 
		<< idx*2+1 << ", " 
		<< idx*12+8 << " , , , "
		<< idx/2+1 << " } 1" << std::endl;
      out << 4*idx+1 << " { " 
		<< idx*2+1 << ", " 
		<< idx*12+6 << " , , , "
		<< "0 } 0" << std::endl;
    }
  }
  
 public:
  tetrahedron_generator(std::ostream& o, RT sin=10) : out(o), s(sin), P(CGAL::to_double(s)/2) {}
    tetrahedron_generator(std::ostream& o, RT sin, unsigned int seed) 
      : out(o), s(sin), r(seed), P(CGAL::to_double(s)/2, r) {}

    void create_tetrahedra(int nx, int ny, int nz) {
      for(int dx=0; dx < nx; ++dx)
	for(int dy=0; dy < ny; ++dy)
	  for(int dz=0; dz < nz; ++dz) 
	    create_tetra(dx,dy,dz);
      print();
    }
};

template <typename K>
const int 
tetrahedron_generator<K>::sprev[6] = {3,5,1,4,0,2};
template <typename K>
const int 
tetrahedron_generator<K>::snext[6] = {4,2,5,0,3,1};
template <typename K>
const int 
tetrahedron_generator<K>::source[6] = {0,1,0,2,1,2};
template <typename K>
const int 
tetrahedron_generator<K>::sface[12] = {0,1,1,0,0,1,
				       1,0,0,1,1,0};
template <typename K>
const int 
tetrahedron_generator<K>::fcycle[8] = {0,1,4,5,3,2,7,6};
template <typename K>
const int 
tetrahedron_generator<K>::facet[24] = {0,1,5,4,2,3,
				    7,6,0,1,5,4,
				    2,3,7,6,0,1,
				    5,4,2,3,7,6};
template <typename K>
const int 
tetrahedron_generator<K>::prev[24] = {8,17,10,19,12,21,
				   14,23,16,1,18,3,
				   20,5,22,7,0,9,
				   2,11,4,13,6,15};
template <typename K>
const int 
tetrahedron_generator<K>::next[24] = {16,9,18,11,20,13,
				   22,15,0,17,2,19,
				   4,21,6,23,8,1,
				   10,3,12,5,14,7};
