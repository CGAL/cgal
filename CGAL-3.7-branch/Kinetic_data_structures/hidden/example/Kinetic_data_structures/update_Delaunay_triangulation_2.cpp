//#define CGAL_CHECK_EXPENSIVE
//#define CGAL_CHECK_EXACTNESS
#define NDEBUG
//#define CGAL_NO_STATIC_FILTERS
//#define MOVE_ALL

#include <CGAL/basic.h>
#include <CGAL/Polynomial/CORE_kernel.h>

#include <CGAL/Random.h>

#include <CGAL/Updatable_Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>

#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <algorithm>


#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

//extern void moncontrol(int);

//#include <cstdio>

template <class TDS>
struct MTDSH
{
  typedef typename TDS::Edge Edge;
  typedef typename TDS::Vertex_handle Vertex_handle;
  typedef typename TDS::Face_handle Face_handle;



  static Edge mirror_edge(const Edge &e) {
    int i= e.first->neighbor(e.second)->index(e.first);
    return Edge(e.first->neighbor(e.second), i);
  }

 
  static Vertex_handle origin(const Edge &e) {
    int o= e.first->ccw(e.second);
    return e.first->vertex(o);
  }

  static Vertex_handle destination(const Edge &e) {
    int o= e.first->cw(e.second);
    return e.first->vertex(o);
  }  

  static Vertex_handle third_vertex(const Edge &e) {
    return e.first->vertex(e.second);
  }
  static Vertex_handle mirror_vertex(const Edge &e) {
    return third_vertex(mirror_edge(e));
  }

  struct Unoriented_edge {
    Unoriented_edge(Face_handle f, int i) {
      if (&*f > &*f->neighbor(i)) {
	e_= mirror_edge(Edge(f,i));
      } else {
	e_ = Edge(f,i);
      }
      
    }
    const Edge& edge() const {
      return e_;
    }
    
    bool operator<(const Unoriented_edge &o) const {
      return &*e_.first < &*o.e_.first || &*e_.first == &*o.e_.first && e_.second < o.e_.second;
    }
    Edge e_;
  };
};



template <class Del>
static bool is_hull_edge(const typename Del::Edge &e) {
  typedef MTDSH <typename Del::Triangulation_data_structure> TDS_helper;
  return ! TDS_helper::mirror_vertex(e)->point().is_valid()
    || ! TDS_helper::third_vertex(e)->point().is_valid()
    || ! TDS_helper::origin(e)->point().is_valid()
    || ! TDS_helper::destination(e)->point().is_valid();
}

template <class Del, class Indirect_kernel>
static bool compute_ok(const typename Del::Edge &e,  Indirect_kernel sk) {
  typedef MTDSH<typename Del::Triangulation_data_structure>  TDS_helper;
  //typename Indirect_kernel::Current_coordinates 
  //cc= sk.current_coordinates_object();
  typename Indirect_kernel::Point_2 ks[4];
  ks[0]= TDS_helper::origin(e)->point();
  ks[1]= TDS_helper::third_vertex(e)->point();
  ks[2]= TDS_helper::destination(e)->point();
  ks[3]= TDS_helper::mirror_vertex(e)->point();
  
  if (ks[0] == typename Indirect_kernel::Point_2() || ks[1]== typename Indirect_kernel::Point_2()
      || ks[2] == typename Indirect_kernel::Point_2() || ks[3] == typename Indirect_kernel::Point_2()){
    typename Indirect_kernel::Orientation_2 o2= sk.orientation_2_object();
  
    
    bool odd_parity=false;
    bool infinity=false;
    for (unsigned int i=0; i<4; ++i) {
      if (infinity) {
	ks[i-1]=ks[i];
      } else {
	if (ks[i] == typename Indirect_kernel::Point_2()) {
	  infinity=true;
	  odd_parity= ((i%2)==1);
	}
      }
    }
    if (odd_parity) {
      std::swap(ks[0], ks[1]);
    }
    CGAL::Orientation o=o2(ks[0], ks[1], ks[2]);
    return o==CGAL::POSITIVE;
  } else {
    typename Indirect_kernel::Side_of_oriented_circle_2 soc
      = sk.side_of_oriented_circle_2_object();
    
    CGAL::Oriented_side s=soc(ks[0], ks[1], ks[2], ks[3]);
    
    if (s== CGAL::ON_ORIENTED_BOUNDARY) {
      std::cout << "Degeneracy with edge " 
		<< ks[0] << " " << ks[2] << std::endl;
    }
    return s!= CGAL::ON_NEGATIVE_SIDE;
  }
}


template <class CNT>
CNT  incircle(CNT ax, CNT ay, CNT bx, CNT by, 
	      CNT cx, CNT cy, CNT dx, CNT dy) {
  CNT qpx = bx - ax;
  CNT qpy = by - ay;
  CNT rpx = cx - ax;
  CNT rpy = cy - ay;
  CNT tpx = dx - ax;
  CNT tpy = dy - ay;
  CNT det=CGAL::determinant(qpx*tpy - qpy*tpx, tpx*(dx - bx) 
				  + tpy*(dy - by),
				  qpx*rpy - qpy*rpx, rpx*(cx - bx) 
				  + rpy*(cy - by));
  return det;
}

template <class Indirect_kernel, class INT>
void ival(typename Indirect_kernel::Point_2 p, typename  Indirect_kernel::Current_coordinates ic,
	  typename Indirect_kernel::Current_coordinates fc, INT rp[2]){
  INT ix= ic(p).x();
  INT iy= ic(p).y();
  INT fx= fc(p).x();
  INT fy= fc(p).y();
  rp[0]= INT(std::min(ix.inf(), fx.inf()), std::max(ix.sup(), fx.sup()));
  rp[1]= INT(std::min(iy.inf(), fy.inf()), std::max(iy.sup(), fy.sup()));
}

template <class Indirect_kernel, class INT>
void ival(typename Indirect_kernel::Point_2 p, 
	  typename Indirect_kernel::Point_2 pb, typename  Indirect_kernel::Current_coordinates ic,
	  typename Indirect_kernel::Current_coordinates fc, INT rp[2]){
  INT ix= ic(p).x();
  INT iy= ic(p).y();
  INT fx= fc(p).x();
  INT fy= fc(p).y();

  INT ixb= ic(pb).x();
  INT iyb= ic(pb).y();
  INT fxb= fc(pb).x();
  INT fyb= fc(pb).y();

  rp[0]= INT(std::min((ix-ixb).inf(), (fx-fxb).inf()), std::max((ix-ixb).sup(), (fx-fxb).sup()));
  rp[1]= INT(std::min((iy-iyb).inf(), (fy-fyb).inf()), std::max((iy-iyb).sup(), (fy-fyb).sup()));
}

template <class Del, class Indirect_kernel>
static bool compute_ok_int(const typename Del::Edge &e,  Indirect_kernel ik, Indirect_kernel fk) {
  typedef MTDSH<typename Del::Triangulation_data_structure>  TDS_helper;
  //typename Indirect_kernel::Current_coordinates 
  //cc= sk.current_coordinates_object();
  typename Indirect_kernel::Point_2 ks[4];
  ks[0]= TDS_helper::origin(e)->point();
  ks[1]= TDS_helper::third_vertex(e)->point();
  ks[2]= TDS_helper::destination(e)->point();
  ks[3]= TDS_helper::mirror_vertex(e)->point();


  
  if (ks[0] == typename Indirect_kernel::Point_2() || ks[1]== typename Indirect_kernel::Point_2()
      || ks[2] == typename Indirect_kernel::Point_2() || ks[3] == typename Indirect_kernel::Point_2()){
    return true;
  } else {
    typename Indirect_kernel::Current_coordinates ic= ik.current_coordinates_object();
    typename Indirect_kernel::Current_coordinates fc= fk.current_coordinates_object();
    typedef CGAL::Interval_nt_advanced INT;
    INT a[2], b[2], c[2], d[2];
    ival<Indirect_kernel>(ks[0], ic, fc, a);
    ival<Indirect_kernel>(ks[1], ic, fc, b);
    ival<Indirect_kernel>(ks[2], ic, fc, c);
    ival<Indirect_kernel>(ks[3], ic, fc, d);
    INT v= incircle(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
    return v.inf() >0;
  }
}



template <class Del, class Indirect_kernel>
static bool compute_ok_intp(const typename Del::Edge &e,  Indirect_kernel ik, Indirect_kernel fk) {
  typedef MTDSH<typename Del::Triangulation_data_structure>  TDS_helper;
  //typename Indirect_kernel::Current_coordinates 
  //cc= sk.current_coordinates_object();
  typename Indirect_kernel::Point_2 ks[4];
  ks[0]= TDS_helper::origin(e)->point();
  ks[1]= TDS_helper::third_vertex(e)->point();
  ks[2]= TDS_helper::destination(e)->point();
  ks[3]= TDS_helper::mirror_vertex(e)->point();


  
  if (ks[0] == typename Indirect_kernel::Point_2() || ks[1]== typename Indirect_kernel::Point_2()
      || ks[2] == typename Indirect_kernel::Point_2() || ks[3] == typename Indirect_kernel::Point_2()){
    return true;
  } else {
    typename Indirect_kernel::Current_coordinates ic= ik.current_coordinates_object();
    typename Indirect_kernel::Current_coordinates fc= fk.current_coordinates_object();
    typedef CGAL::Interval_nt_advanced INT;
    INT b[2], c[2], d[2];
    //ival<Indirect_kernel>(ks[0], ic, fc, a);
    ival<Indirect_kernel>(ks[1], ks[0], ic, fc, b);
    ival<Indirect_kernel>(ks[2], ks[0], ic, fc, c);
    ival<Indirect_kernel>(ks[3], ks[0], ic, fc, d);
    INT v= incircle(INT(0.0), INT(0.0), b[0], b[1], c[0], c[1], d[0], d[1]);
    return v.inf() >0;
  }
}




template <class Points>
void read(std::string name, Points &points) {
    std::ifstream in(name.c_str());
    if (!in) {
      std::cerr << "Error opening file " << name  << std::endl;
      exit(1);
    }
    
    while (true) {
      char ln[10000];
      in.getline(ln, 10000);
      if (!in) {
	break;
      }
      typename Points::value_type p;
      std::istringstream iss(ln);
      iss >> p;
      if (!iss) {
	CGAL::Simple_cartesian<double>::Point_2 dpt;
	std::istringstream iss2(ln);
	iss2 >> dpt;
	if (!iss2) {
	  std::cerr << "Error processing line " << ln << std::endl;
	} else {
	  points.push_back(typename Points::value_type(dpt.x(), dpt.y()));
	}
      } else {

	points.push_back(p);
      }
    };
  }

template <class Tri, class IK>
int hybrid(Tri &del, IK, IK fk) {
  typedef MTDSH<typename Tri::Triangulation_data_structure>  TDSH;
  typedef typename Tri::Edge Edge;
  
  int bc=0;
  typename IK::Orientation_2 o2= fk.orientation_2_object();
  for (typename Tri::Finite_faces_iterator fit = del.finite_faces_begin(); fit != del.finite_faces_end(); ++fit){
    if (o2(fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point()) != CGAL::POSITIVE) {
      ++bc;
    }
  }
  if (bc != 0) {
    std::cout << "bc is not zero " << std::endl;
  }
  typename std::set<typename TDSH::Unoriented_edge> bad;
  //typename IK::Side_of_oriented_circle_2 c2= fk.side_of_oriented_circle_2_object();
  for (typename Tri::Finite_edges_iterator fit = del.finite_edges_begin(); fit != del.finite_edges_end(); ++fit){
    if (!compute_ok<Tri, IK>(*fit, fk)) {
      bad.insert(typename TDSH::Unoriented_edge(fit->first, fit->second));
    }
  }

  int nflips= 0;
  while (!bad.empty()) {
    Edge e= bad.begin()->edge();
    bad.erase(bad.begin());
    Edge me= TDSH::mirror_edge(e);
    ++nflips;
    for (int i=0; i< 3; ++i) {
      if (i != e.second) {
	bad.erase(typename TDSH::Unoriented_edge(e.first, i));
      }
      if (i != me.second) {
	bad.erase(typename TDSH::Unoriented_edge(me.first, i));
      }
    }
   
    del.flip(e.first, e.second);
  
    for (int i=0; i< 3; ++i) {
      if (i != e.second) {
	if (!compute_ok<Tri, IK>(Edge(e.first, i), fk)) {
	  bad.insert(typename TDSH::Unoriented_edge(e.first, i));
	}
      }
      if (i != me.second) {
	if (!compute_ok<Tri, IK>(Edge(me.first, i), fk)) {
	  bad.insert(typename TDSH::Unoriented_edge(me.first, i));
	}
      }
    }
  }
  //del.geom_traits().swap(fk);
  return nflips;
}

int main(int argc, char *argv[]) {
  
  //moncontrol(0);

  bool print_help=false;
  bool speedup_only=false;
  bool profile=true;
  std::string ifile="../../demo/Kinetic_data_structures/data/before106", 
    ffile="../../demo/Kinetic_data_structures/data/after106";
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("speedup-only,s", boost::program_options::bool_switch(&speedup_only), "Only show the speedup and the fraction failed")
    ("disable-0,0", boost::program_options::bool_switch(&CGAL::disable_filter_0_), "Disable filter 0")
    ("disable-1,1", boost::program_options::bool_switch(&CGAL::disable_filter_1_), "Disable filter 1")
    ("disable-2,2", boost::program_options::bool_switch(&CGAL::disable_filter_2_), "Disable filter 2")
    ("disable-3,3", boost::program_options::bool_switch(&CGAL::disable_filter_3_), "Disable filter 3")

    ("profile,p", boost::program_options::bool_switch(&profile), "Just update for profiling")

    ("ifile,i", boost::program_options::value<std::string>(&ifile), "The inital coordinates.")
    ("ffile,f", boost::program_options::value<std::string>(&ffile), "The final coordinates.")
    ;

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(desc).run(), vm);
  boost::program_options::notify(vm);

  if (print_help) {
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }
#endif
  
  //if (exact) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel  CK;
  typedef CGAL::Indirect_point_2_kernel<CGAL::Exact_predicates_inexact_constructions_kernel>  K;
  typedef CGAL::Updatable_Delaunay_triangulation_2<K> UD;
  typedef CGAL::Delaunay_triangulation_2<K> Del;
  typedef CGAL::Delaunay_triangulation_2<CK> CDel;


  typedef CGAL::Triangulation_vertex_base_2<K>             HVbb;
  typedef CGAL::Triangulation_hierarchy_vertex_base_2<HVbb> HVb;
  typedef CGAL::Triangulation_face_base_2<K>               HFb;
  typedef CGAL::Triangulation_data_structure_2<HVb,HFb>      HTds;
  typedef CGAL::Delaunay_triangulation_2<K,HTds>            HDt;
  typedef CGAL::Triangulation_hierarchy_2<HDt>              H;
  
  typedef UD::Points Points;
  Points ipts, fpts;
  read(ifile, ipts);
  read(ffile, fpts);
  
  {
    K ik, fk;
    K::Key_range ir= ik.new_point_2s(ipts.begin(), ipts.end());
    K::Key_range fr= fk.new_point_2s(fpts.begin(), fpts.end());
    
    Del d(ik);
    d.insert(ir.first, ir.second);
    int bad_count=0;
    int filter_1_failure=0;
    int filter_1p_failure=0;
    for (Del::Edge_iterator it= d.edges_begin(); 
	 it != d.edges_end(); ++it){
      if (!compute_ok<Del>(*it, fk)) {
	++bad_count;
      }
      if (!compute_ok_int<Del>(*it, ik, fk)) {
	++filter_1_failure;
      }
      if (!compute_ok_intp<Del>(*it, ik, fk)) {
	++filter_1p_failure;
      }
    }
    if (!speedup_only) {
      std::cout << "For " << std::distance(d.edges_begin(), d.edges_end()) << " edges " 
		<< " there were " << bad_count << " failed certs "
		<< "and " <<  filter_1_failure << " edges which failed on 1" 
		<< " and " << filter_1p_failure << " which failed on the improved 1" << std::endl;
    } else {
      std::cout << static_cast<double>(bad_count)/std::distance(d.edges_begin(), d.edges_end()) << " ";
    }
  }
  int update_speed;
  int activated=0;
  {
    UD ud(ipts.begin(), ipts.end());

    // moncontrol(1);

    if (!speedup_only) std::cout << "Timing patching..." << std::flush;
    
    int count=0;
    double min_time=10;
    if (profile) min_time =60;
    CGAL::Timer tim;
    tim.start();
    while (tim.time() < min_time) {
      ud.update_coordinates(fpts);
      if (count ==0 ) activated= CGAL::stat_number_of_activated_vertices_;
      ud.update_coordinates(ipts);
      count+=2;
    }
    tim.stop();
    if (!speedup_only) {
      std::cout << "done." << std::endl;
      std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
      ud.write_statistics(std::cout);
    }
    ud.audit();
    update_speed= count;
  }
  
  if (!profile) {
    K ik, fk;
    K::Key_range ir= ik.new_point_2s(ipts.begin(), ipts.end());
    fk.new_point_2s(fpts.begin(), fpts.end());
    if (!speedup_only) std::cout << "Timing hybrid..." << std::flush;
    Del d(ik);
    d.insert(ir.first, ir.second);
    int count=0; 
    CGAL::Timer tim;
    tim.start();
    int nflips;
    do  {
      nflips= hybrid(d, ik, fk);
      nflips+= hybrid(d, fk, ik);
      //std::cout << nflips << std::endl;
      count+=2;
   } while (tim.time() < 10);
    tim.stop();
    if (!speedup_only) {
      std::cout << "done." << std::endl;
      std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
      std::cout << "nflips is " << nflips << std::endl;
    } 
  }
  if (!profile) {
    K ik, fk;
    K::Key_range ir= ik.new_point_2s(ipts.begin(), ipts.end());
    K::Key_range fr= fk.new_point_2s(fpts.begin(), fpts.end());
    if (!speedup_only) std::cout << "Timing verifying..." << std::flush;
    
    Del d(ik);
    d.insert(ir.first, ir.second);

    int count=0; 
    CGAL::Timer tim;
    tim.start();
    int bad_count;
    do  {
      bad_count=0;
      for (Del::Edge_iterator it= d.edges_begin(); 
	   it != d.edges_end(); ++it){
	if (!compute_ok<Del>(*it, fk)) {
	  ++bad_count;
	}
      }
      ++count;
   } while (tim.time() < 10);
    tim.stop();
    if (!speedup_only) {
      std::cout << "done." << std::endl;
      std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
      std::cout << "Bad count is " << bad_count << std::endl;
    } 
  }
  
  if (0) {
    K ik, fk;
    K::Key_range ir= ik.new_point_2s(ipts.begin(), ipts.end());
    K::Key_range fr= fk.new_point_2s(fpts.begin(), fpts.end());
    std::cout << "Timing rebuilding..." << std::flush;
    std::vector<K::Point_2> ipts(ir.first, ir.second);

    int count=0; 
    CGAL::Timer tim;
    tim.start();
    while (tim.time() < 10) {
      std::random_shuffle(ipts.begin(), ipts.end());
      Del d(ik);
      d.insert(ipts.begin(), ipts.end());
      std::random_shuffle(ipts.begin(), ipts.end());
      Del d2(fk);
      d2.insert(ipts.begin(), ipts.end());
      count+=2;
    }
    tim.stop();
    std::cout << "done." << std::endl;
    std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
  }

  if (!profile) {
    K ik, fk;
    K::Key_range ir= ik.new_point_2s(ipts.begin(), ipts.end());
    K::Key_range fr= fk.new_point_2s(fpts.begin(), fpts.end());
    if (!speedup_only) std::cout << "Timing rebuilding with hierarchy..." << std::flush;
    std::vector<K::Point_2> ipts(ir.first, ir.second);

    int count=0; 
    CGAL::Timer tim;
    tim.start();
    while (tim.time() < 10) {
      std::random_shuffle(ipts.begin(), ipts.end());
      H d(ik);
      d.insert(ipts.begin(), ipts.end());
      std::random_shuffle(ipts.begin(), ipts.end());
      H d2(fk);
      d2.insert(ipts.begin(), ipts.end());
      count+=2;
    }
    tim.stop();
    if (!speedup_only) {
      std::cout << "done." << std::endl;
      std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
    } else {
      std::cout << static_cast<double>(update_speed)/count << " " << static_cast<double>(activated)/
	CGAL::stat_number_of_vertices_ << std::endl;
    }
  }

  if (!profile) {
    CDel d;
    std::vector<CDel::Vertex_handle> vhs(ipts.size());
    for (unsigned int i=0; i < ipts.size(); ++i) {
      vhs[i]= d.insert(ipts[i]);
      if (vhs[i] == CDel::Vertex_handle()) {
	std::cerr  << "Point " << ipts[i] << " duplicated in input." << std::endl;
      }
    }
    if (0) {
      std::cout << "Timing reinsert ..." << std::flush;
      
      int count=0; 
      CGAL::Timer tim;
      
      tim.start();
      while (tim.time() < 10) {
	for (unsigned int i=0; i < ipts.size(); ++i) {
	  if (fpts[i] != vhs[i]->point()) {
	    CDel::Vertex_handle nvhs= d.insert(fpts[i], vhs[i]->face());
	    d.remove(vhs[i]);
	    vhs[i]=nvhs;
	  }
	}
	for (unsigned int i=0; i < fpts.size(); ++i) {
	  if (ipts[i] != vhs[i]->point()) {
	    CDel::Vertex_handle nvhs=  d.insert(ipts[i], vhs[i]->face());
	    d.remove(vhs[i]);
	    vhs[i]=nvhs;
	  }
	}
	count+=2;
      }
      tim.stop();
      std::cout << "done." << std::endl;
      std::cout << count << " reps in " << tim.time() 
		<< " seconds" << std::endl;
    }
  }
  if (0){
    UD ud(ipts.begin(), ipts.end());

    // moncontrol(1);

    std::cout << "Timing patching..." << std::flush;

    int count=0;
    double min_time=10;
    if (profile) min_time =60;
    CGAL::Timer tim;
    tim.start();
    while (tim.time() < min_time) {
      ud.update_coordinates(fpts);
      ud.update_coordinates(ipts);
      count+=2;
    }
    tim.stop();

    std::cout << "done." << std::endl;
    std::cout << count << " reps in " << tim.time() << " seconds" << std::endl;
    ud.write_statistics(std::cout);
    ud.audit();
  }


  return EXIT_SUCCESS;
  
};
