//#define  CGAL_DISABLE_AM_CODE 1
//#define  CGAL_DISABLE_M_CODE 1
// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// add this to profile
//#define CGAL_PROFILE 1

// add this to dump quadruples of sites for which the incircle test is called
//#define CGAL_PROFILE_SDG_DUMP_INCIRCLE 1

// add this to test if both old and new incircle tests give the same result
//#define CGAL_SDG_CHECK_INCIRCLE_CONSISTENCY 1

// add this to use the old incircle test
//#define CGAL_SDG_USE_OLD_INCIRCLE 1

#define CGAL_SDG_NO_FACE_MAP 1
#define USE_INPLACE_LIST 1
#define CGAL_SDG_USE_SIMPLIFIED_ARRANGEMENT_TYPE_PREDICATE 1
//#define CGAL_SDG_SORT_POINTS_IN_SITE2 1

// choose the kernel
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/Leda_real.h>
#else
//#  include <CGAL/CORE_Expr.h>
#  include <CGAL/Gmpq.h>
#endif


typedef CGAL::Simple_cartesian<double> K;
typedef  K::Point_2 Point_2;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>

typedef CGAL::Field_with_sqrt_tag MTag;
#ifdef CGAL_USE_LEDA
typedef CGAL::Field_with_sqrt_tag EMTag;
typedef CGAL::Simple_cartesian<leda::real> EK;
#else
//typedef CGAL::Field_with_sqrt_tag EMTag;
//typedef CGAL::Simple_cartesian<CORE::Expr> EK;
typedef CGAL::Integral_domain_without_division_tag EMTag;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> EK;
#endif
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<K, MTag, EK, EMTag>  Gt;


#ifdef USE_INPLACE_LIST

typedef CGAL::Segment_Delaunay_graph_storage_traits_2<Gt> ST;
typedef CGAL::Segment_Delaunay_graph_vertex_base_2<ST>    Vb;
typedef CGAL::Segment_Delaunay_graph_face_base_2<Gt>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Segment_Delaunay_graph_2<Gt,ST,TDS,CGAL::Tag_true>  SDG2;
#else
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG2;
#endif

typedef std::vector<Point_2> Points_container;

typedef  Points_container::size_type Index_type;
typedef std::pair<Index_type,Index_type> Constraint;

typedef std::vector<Constraint> Constraints_container;
Points_container points;
Constraints_container constraints;
SDG2 sdg;

std::vector<SDG2::Site_2> sites;



template <typename Kernel, typename Iterator>
struct Sort_traits_2 {
  Kernel k;
  Sort_traits_2 (const Kernel &kernel = Kernel())
    : k (kernel)
  {}

  typedef Iterator Point_2;

  struct Less_x_2 {
    Kernel k;
    Less_x_2 (const Kernel &kernel = Kernel())
      : k (kernel)
    {}
    bool operator() (const Point_2 &p, const Point_2 &q) const
    {
      return k.less_x_2_object() (*p, *q);
    }
  };

  Less_x_2
  less_x_2_object() const
  {
    return Less_x_2(k);
  }

  struct Less_y_2 {
    Kernel k;
    Less_y_2 (const Kernel &kernel = Kernel())
      : k (kernel)
    {}
    bool operator() (const Point_2 &p, const Point_2 &q) const
    {
      return k.less_y_2_object() (*p, *q);
    }
  };

  Less_y_2
  less_y_2_object() const
  {
    return Less_y_2(k);
  }
};


template <typename SDG>
void
insert_constraints_using_spatial_sort(SDG& sdg)
{
  typedef typename Points_container::const_iterator Points_iterator;
  typedef std::vector<Points_iterator> Indices;
  typedef std::vector<typename SDG::Vertex_handle> Vertices;

  Sort_traits_2<K, Points_iterator> sort_traits;

  Indices indices;
  indices.reserve(points.size());
  for(Points_iterator it = points.begin(); it != points.end(); ++it) {
    indices.push_back(it);
  }
  std::cerr << "Spatial sorting...";
  CGAL::Timer timer;
  timer.start();
  std::random_shuffle(indices.begin(), indices.end());
  CGAL::spatial_sort(indices.begin(), indices.end(),
                     sort_traits);
  timer.stop();
  std::cerr << " done (" << timer.time() << "s)\n";

  std::cerr << "Inserting points...";
  timer.reset();
  timer.start();
  Vertices vertices;
  vertices.resize(points.size());
  typename SDG::Vertex_handle hint;
  for(typename Indices::const_iterator
        pt_it_it = indices.begin(), end = indices.end();
      pt_it_it != end; ++pt_it_it) {
    typename SDG::Vertex_handle vh = sdg.insert(**pt_it_it, hint);
    hint = vh;
    vertices[*pt_it_it - points.begin()] = vh;
  }
  timer.stop();
  std::cerr << " done (" << timer.time() << "s)\n";

  std::cerr << "Inserting constraints...";

  timer.reset();
  timer.start();
  for(typename Constraints_container::const_iterator
        cit = constraints.begin(), end = constraints.end();
      cit != end; ++cit) {
    const typename SDG::Vertex_handle& v1 = vertices[cit->first];
    const typename SDG::Vertex_handle& v2 = vertices[cit->second];
    if(v1 != v2)
      sdg.insert(v1, v2);
  }

  timer.stop();
  std::cerr << " done (" << timer.time() << "s)\n";
}


bool
load_cin_file(std::istream& ifs) {
  std::cerr << "Loading file... ";
  CGAL::Timer timer;
  timer.start();
  bool not_first = false;
  int n;
  ifs >> n;
  if(!ifs.good())
    return false;
  Point_2 p, q, qold;
  char c;
  int point_counter = 0;
  while(ifs >> c >> p >> q) {
    if(not_first && p == q) {
      continue;
    }
    if(p == qold) {
      points.push_back(q);
      constraints.push_back(std::make_pair(point_counter-1, point_counter));
      ++point_counter;
    }
    else {
      points.push_back(p);
      points.push_back(q);
      constraints.push_back(std::make_pair(point_counter, point_counter+1));
      point_counter += 2;
    }
    qold = q;
    not_first = true;
  }
  std::cerr << "done (" << timer.time() << "s)" << std::endl;
  insert_constraints_using_spatial_sort(sdg);
  return true;
}


int main()
{
  load_cin_file(std::cin);
  if(! sdg.is_valid(true, 1) ){
    std::cerr << "invalid data structure" << std::endl;
  } else {
    std::cerr << "valid data structure" << std::endl;
  }

  return 0;
}

