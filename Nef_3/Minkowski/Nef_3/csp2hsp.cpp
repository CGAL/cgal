#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#else
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_3/Bounding_box_3.h>
#include <CGAL/Nef_3/Nary_union.h>
#include <CGAL/Nef_3/volInt.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_3/convex_decomposition_3.h>
#include <CGAL/convexity_check_3.h>

#include <fstream>
#include <sstream>

//#define CGAL_WITH_LAZY_KERNEL
#ifdef CGAL_WITH_LAZY_KERNEL
typedef CGAL::Gmpq NT;
//typedef leda_rational NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
#else
#ifdef CGAL_USE_LEDA
typedef leda_integer NT;
#else
typedef CGAL::Gmpz NT;
#endif
typedef CGAL::Homogeneous<NT> Kernel;
#endif
typedef Kernel::RT RT;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
//typedef Polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
#ifdef CGAL_NEF_INDEXED_ITEMS
typedef CGAL::Nef_polyhedron_3<Kernel,CGAL::SNC_indexed_items>     Nef_polyhedron_3;
#else
typedef CGAL::Nef_polyhedron_3<Kernel>     Nef_polyhedron_3;
#endif
typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef Nef_polyhedron_3::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron_3::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron_3::Halffacet_const_iterator Halffacet_const_iterator;
typedef Nef_polyhedron_3::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron_3::Volume_const_handle Volume_const_handle;
typedef Nef_polyhedron_3::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron_3::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron_3::SFace_const_handle SFace_const_handle;
typedef CGAL::Bounding_box_3<CGAL::Tag_true, Kernel> BBox;
typedef CGAL::Nary_union<Nef_polyhedron_3> Nary_union;

class Volume_output {
  bool twin;
  std::vector<Plane_3> planes;
public:

  typedef std::vector<Plane_3>::const_iterator CI;

  Volume_output(bool twin_ = true) : twin(twin_){}
  bool is_in(const Plane_3 p) {
    for(CI ci = planes.begin(); ci != planes.end(); ++ci)
      if(*ci == p)
	return true;
    return false;
  }

  void visit(Halffacet_const_handle f) {
    if(twin) f = f->twin();
    Plane_3 p = f->twin()->plane();
    if(!is_in(p))
      planes.push_back(p);
  }
  void visit(SFace_const_handle s) {}
  void visit(Halfedge_const_handle e) {}
  void visit(Vertex_const_handle v) {}
  void visit(SHalfedge_const_handle se) {}
  void visit(SHalfloop_const_handle sl) {}

  void dump() const {
    std::cout << planes.size() << std::endl;
    for(CI ci = planes.begin(); ci != planes.end(); ++ci) {
      Plane_3 p = *ci;

#ifdef CGAL_WITH_LAZY_KERNEL
      std::cout << CGAL::to_double(p.a().exact()) << " "
                << CGAL::to_double(p.b().exact()) << " "
                << CGAL::to_double(p.c().exact()) << " "
                << CGAL::to_double(p.d().exact()) << std::endl;
#else
#ifdef CGAL_CSP2HSP_EXACT_OUTPUT
      std::cout << p << std::endl;
#else
      std::cout << CGAL::to_double(p.a()) << " "
                << CGAL::to_double(p.b()) << " "
                << CGAL::to_double(p.c()) << " "
                << CGAL::to_double(p.d()) << std::endl;
#endif
#endif
    }
  }
    
  CI planes_begin() const {
    return planes.begin();
  }
    
  CI planes_end() const {
    return planes.end();
  }
};


class Point_collector {

  std::list<Point_3> points;
public:

  Point_collector() {}

  void visit(Vertex_const_handle v) {
    points.push_back(v->point());
  }

  void visit(Halffacet_const_handle) {}
  void visit(SFace_const_handle) {}
  void visit(Halfedge_const_handle) {}
  void visit(SHalfedge_const_handle) {}
  void visit(SHalfloop_const_handle) {}

  std::list<Point_3>::const_iterator points_begin() const {
    return points.begin();
  }
    
  std::list<Point_3>::const_iterator points_end() const {
    return points.end();
  }
};


class Convex_from_shell_points {

  std::list<Point_3> points;

public:
  Convex_from_shell_points() {}

  void visit(Vertex_const_handle v) {
    points.push_back(v->point());
  }

  void visit(Halffacet_const_handle) {}
  void visit(SFace_const_handle) {}
  void visit(Halfedge_const_handle) {}
  void visit(SHalfedge_const_handle) {}
  void visit(SHalfloop_const_handle) {}

  Nef_polyhedron_3 get_polyhedron() 
  { 
    Polyhedron_3 cv;
    convex_hull_3(points.begin(), points.end(), cv);
    return Nef_polyhedron_3(cv); 
  }

  void clear() { points.clear(); }

};

class BBox_constructor {

  BBox bbox;
public:

  BBox_constructor() : bbox() {}

  void visit(Vertex_const_handle v) {
    bbox.extend(v->point());
  }

  void visit(Halffacet_const_handle) {}
  void visit(SFace_const_handle) {}
  void visit(Halfedge_const_handle) {}
  void visit(SHalfedge_const_handle) {}
  void visit(SHalfloop_const_handle) {}

  BBox get_box() { return bbox; }
};


template <typename CI>
Nef_polyhedron_3 create_from_halfspaces(CI begin, CI end, bool invert) {
  
  Nef_polyhedron_3 cube;
  std::ifstream in("Nef3/centered_cube.nef3");
  in >> cube;

  int mag = 1000000;
  Aff_transformation_3 scale(mag, 0, 0, 
			     0, mag, 0, 
			     0, 0, mag, 1);
  cube.transform(scale);
 
  CI pli;
  for(pli = begin; pli != end; ++pli) {
    if(invert) {
      Plane_3 inverted_plane(-pli->a(), -pli->b(),
			     -pli->c(), -pli->d());
      cube = cube.intersection(inverted_plane);
    } else
      cube = cube.intersection(*pli);
  }  
  return cube;
}

void output_hsp(const Nef_polyhedron_3& DIFF,
		const Nef_polyhedron_3& NCV) {

  int nov = 1;
  Volume_const_iterator c;
  for(c=++(DIFF.volumes_begin()); c!=DIFF.volumes_end(); ++c)
    if(c->mark()) ++nov;
  std::cout << nov << std::endl;

  Volume_output vout_ncv;
  NCV.visit_shell_objects(NCV.volumes_begin()->shells_begin(), vout_ncv);
  vout_ncv.dump();

  for(c=++(DIFF.volumes_begin()); c!=DIFF.volumes_end(); ++c) {
    if(c->mark()) {
      Volume_output vout(true);
      DIFF.visit_shell_objects(c->shells_begin(), vout);
      vout.dump();     
 
#ifdef CGAL_CSP2HSP_TEST_CONVEX_PARTS
      Nef_polyhedron_3 obs1 = 
	create_from_halfspaces(vout.planes_begin(), vout.planes_end(), true);

      Polyhedron_3 pobs2;
      DIFF.convert_inner_shell_to_polyhedron(c->shells_begin(), pobs2);
      Nef_polyhedron_3 obs2(pobs2);      

      Nef_polyhedron_3 empty = obs1.symmetric_difference(obs2);
      if(!empty.is_empty())
	std::cerr << "error " << std::endl;
#endif
    }
  }
}

/*
double bbox_size(const BBox& bbox) {

  FT size =
    (bbox.max_coord(0) - bbox.min_coord(0))*
    (bbox.max_coord(1) - bbox.min_coord(1))*
    (bbox.max_coord(2) - bbox.min_coord(2));
  return CGAL::to_double(size);
}
*/

template<typename Hash_map>
Volume_const_handle find(Volume_const_handle c, Hash_map& hash) {
  if(c != hash[c])
    hash[c] = find(hash[c], hash);
  return hash[c];
}

void simplify_and_output(const Nef_polyhedron_3& DIFF,
			 const Nef_polyhedron_3& NCV,
			 int simplify, double factor, int simplifyBy, double minsize) {
  
  typedef CGAL::PolyhedralVolumeCalculator<Nef_polyhedron_3> PVC;
  typedef CGAL::Union_find<Volume_const_handle> UF;
  typedef UF::handle UF_handle;
  typedef UF::iterator UF_iterator;
  typedef UF::const_pointer UF_pointer;
  UF uf;
  
  CGAL::Unique_hash_map<Volume_const_handle, Volume_const_handle> c2c;
  CGAL::Unique_hash_map<Volume_const_handle, Nef_polyhedron_3> c2N;
  CGAL::Unique_hash_map<Volume_const_handle, double> real_volume;

  PVC pvc(DIFF);

  int volumes = 0;
  Volume_const_iterator ci;
  CGAL_forall_volumes(ci, DIFF) {
    if(!ci->mark()) continue;
    ++volumes;
    c2c[ci] = ci;
    Convex_from_shell_points cfsp;
    DIFF.visit_shell_objects(ci->shells_begin(), cfsp);
    c2N[ci] = cfsp.get_polyhedron();
    real_volume[ci] = pvc.get_volume_of_closed_volume(ci);
  }

  bool simplified;
  do{
    simplified = false;
    std::cerr << "next run " << std::endl;

    Halffacet_const_iterator fi;
    CGAL_forall_halffacets(fi, DIFF) {
      if(fi->is_twin()) continue;
      Volume_const_handle c0 = fi->incident_volume();
      Volume_const_handle c1 = fi->twin()->incident_volume();
      if(!c0->mark() || !c1->mark()) continue;
      c0 = find(c0, c2c);
      c1 = find(c1, c2c);
      if(c0 == c1) continue;
      
      Nef_polyhedron_3 N0 = c2N[c0];
      Nef_polyhedron_3 N1 = c2N[c1];
      
      std::list<Point_3> points;
      Vertex_const_iterator vi;
      CGAL_forall_vertices(vi, N0)
	points.push_back(vi->point());
      CGAL_forall_vertices(vi, N1)
	points.push_back(vi->point());
      
      Polyhedron_3 cv;
      convex_hull_3(points.begin(), 
		    points.end(), cv);
      Nef_polyhedron_3 ncv(cv);
      
      PVC pvct(ncv);
      double s0 = real_volume[c0];
      double s1 = real_volume[c1];
      double s01 = pvct.get_volume_of_polyhedron();
      CGAL_assertion(s0+s1 <= s01);
      if(simplify < 1 || simplify > 3) continue;
      if(((simplify & 1) == 0 || (s0+s1)*factor < s01) &&
	 ((simplify & 2) == 0 || s0+s1+simplifyBy < s01)) continue;
      
      c2c[c1] = c0;
      c2N[c0] = ncv;
      real_volume[c0] = s0+s1;
      --volumes;
      simplified = true;
    }

    std::cerr << "obstacles left " << volumes << std::endl;
  } while(simplified);

  Nary_union nu2;
  CGAL_forall_volumes(ci, DIFF) {
    if(find(ci, c2c) != ci) continue;
    nu2.add_polyhedron(c2N[ci]);
  }

  Nef_polyhedron_3 rediff = nu2.get_union();
  Nef_polyhedron_3 recv2 = NCV - rediff;  
  PVC pvct2(recv2);
  std::cerr << "volume of approximation " 
	    << pvct2.get_volume_of_polyhedron() << std::endl;
  CGAL_assertion((recv2-NCV).is_empty());

  if(minsize > 0.0) {
    std::cerr << "kickout obstacles smaller than " << minsize << "cmm" << std::endl;
    CGAL_forall_volumes(ci, DIFF) {
      if(find(ci, c2c) != ci) continue;
      if(minsize<=0 || real_volume[ci] < minsize) --volumes; 
    }
    std::cerr << "obstacles left " << volumes << std::endl;
  }

  std::cout << volumes+1 << std::endl;

  Volume_output vout_ncv;
  NCV.visit_shell_objects(NCV.volumes_begin()->shells_begin(), vout_ncv);
  vout_ncv.dump();

  CGAL_forall_volumes(ci, DIFF) {
    if(find(ci, c2c) != ci) continue;
    if(minsize>0 && real_volume[ci] < minsize) continue; 
    Nef_polyhedron_3 ncv = c2N[ci];
    Volume_output vout;
    ncv.visit_shell_objects(ncv.volumes_begin()->shells_begin(), vout);
    vout.dump();
  }

}

int main(int argc, char* argv[]) {

  int simplify = argc > 1 ? std::atoi(argv[1]) : 0;
  int prozent = argc > 2 ? std::atoi(argv[2]) : 20;
  int simplifyBy = argc > 3 ? std::atoi(argv[3]) : 30000;
  double minsize = argc > 4 ? std::atof(argv[4]) : 0.0;
    
  double factor = (100.0+double(prozent))/100.0;

  if(simplify > 0)
    std::cerr << "simplify " << simplify << ", " 
	      << factor << ", " << simplifyBy << std::endl;

  Nef_polyhedron_3 CSP;
  std::cin >> CSP;

  if(CSP.number_of_vertices()==0)
    return 0;

  std::list<Point_3> points;
#ifdef CGAL_CSP2HSP_BOX_AS_CV
  Vertex_const_iterator v = CSP.vertices_begin();
  FT p[3];
  p[0] = v->point().x();
  p[1] = v->point().y();
  p[2] = v->point().z();
  BBox B(p);
  for(++v; v != CSP.vertices_end(); ++v)
    B.extend(v->point());

  points.push_back(Point_3(B.min_coord(0), 
			B.min_coord(1),
			B.min_coord(2)));
  points.push_back(Point_3(B.min_coord(0), 
			B.min_coord(1),
			B.max_coord(2)));
  points.push_back(Point_3(B.min_coord(0), 
			B.max_coord(1),
			B.min_coord(2)));
  points.push_back(Point_3(B.min_coord(0), 
			B.max_coord(1),
			B.max_coord(2)));
  points.push_back(Point_3(B.max_coord(0), 
			B.min_coord(1),
			B.min_coord(2)));
  points.push_back(Point_3(B.max_coord(0), 
			B.min_coord(1),
			B.max_coord(2)));
  points.push_back(Point_3(B.max_coord(0),
			B.max_coord(1),
			B.min_coord(2)));
  points.push_back(Point_3(B.max_coord(0), 
			B.max_coord(1),
			B.max_coord(2))); 

#else
  Vertex_const_iterator v;
  for(v = CSP.vertices_begin(); v != CSP.vertices_end(); ++v)
    points.push_back(v->point());
#endif

  Polyhedron_3 CV;
  convex_hull_3( points.begin(), points.end(), CV);

  Nef_polyhedron_3 NCV(CV);
  Nef_polyhedron_3 DIFF = NCV-CSP.interior();

  typedef CGAL::PolyhedralVolumeCalculator<Nef_polyhedron_3> PVC;
  PVC pvc_ncv(NCV);
  double sncv = pvc_ncv.get_volume_of_polyhedron();
  PVC pvc_diff(DIFF);
  double sdiff = pvc_diff.get_volume_of_polyhedron();
  std::cerr << "ncv - diff = " << sncv << " - " << sdiff 
	    << " = " << sncv-sdiff << std::endl;
  
  convex_decomposition_3<Nef_polyhedron_3>(DIFF);

  std::ofstream outDiff("diff.nef3");
  outDiff << DIFF;

  if(simplify > 0)
    simplify_and_output(DIFF, NCV, simplify, factor, simplifyBy, minsize);
  else
    output_hsp(DIFF, NCV);
}
