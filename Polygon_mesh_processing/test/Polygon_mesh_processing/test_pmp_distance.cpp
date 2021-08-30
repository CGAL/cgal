#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Real_timer.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/number_utils.h>
#include <CGAL/Coercion_traits.h>

#include <fstream>
#include <ostream>

struct Custom_traits_Hausdorff
{
// New requirements {
  struct FT
  {
    FT(){}
    FT(double){}
    FT operator/(FT)const{return FT();}
    FT operator-(const FT)const{return FT();}
    FT operator+(FT)const{return FT();}
    FT operator*(FT)const{return FT();}
    bool operator<=(FT)const{return false;}
    bool operator>=(FT)const{return false;}
    bool operator!=(FT)const{return false;}
    bool operator<(FT)const{return false;}
    bool operator>(FT)const{return false;}
    FT& operator+=(FT){ return *this; }
  };

  struct Point_3
  {
    Point_3(){}
    Point_3(FT, FT, FT){}
    FT operator[](int)const{return FT();}
    bool operator<(Point_3)const{return false;}
    double x()const{return 0;}
    double y()const{return 0;}
    double z()const{return 0;}
  };

  struct Vector_3{};

  struct Triangle_3
  {
    Triangle_3(const Point_3&,  const Point_3&, const Point_3&){}
    Point_3 operator[](int)const{return Point_3();}
    CGAL::Bbox_3 bbox(){return CGAL::Bbox_3();}
  };

  struct Segment_3
  {
    Segment_3(const Point_3&,  const Point_3&){}
    Point_3 operator[](int)const{return Point_3();}
  };

  struct Compute_squared_area_3
  {
    typedef FT result_type;
    FT operator()(Point_3,Point_3,Point_3)const{return FT();}
    FT operator()(Triangle_3)const{return FT();}
  };

  struct Compute_squared_length_3
  {
    typedef FT result_type;
    FT operator()(Segment_3)const{return FT();}
  };

  struct Construct_translated_point_3
  {
    Point_3 operator() (const Point_3 &, const Vector_3 &){return Point_3();}
  };

  struct Construct_vector_3{
    Vector_3         operator() (const Point_3 &, const Point_3 &){return Vector_3();}
  };

  struct Construct_scaled_vector_3
  {
    Vector_3         operator() (const Vector_3 &, const FT &)
    {return Vector_3();}

  };

  Compute_squared_area_3 compute_squared_area_3_object(){return Compute_squared_area_3();}
// } end of new requirements

// requirements from AABBGeomTraits {
  struct Sphere_3
  {};

  struct Equal_3
  {
    bool operator()(Point_3,Point_3) const {return false;}
  };

  struct Do_intersect_3
  {
    CGAL::Comparison_result operator()(const Sphere_3& , const CGAL::Bbox_3&){ return CGAL::ZERO;}
  };

  struct Intersect_3
  {
    struct result_type{};
  };

  struct Construct_sphere_3
  {
    Sphere_3 operator()(const Point_3& , FT ){return Sphere_3();}
  };

  struct Construct_projected_point_3
  {
    const Point_3 operator()(Triangle_3,  Point_3)const{return Point_3();}
  };

  struct Compare_distance_3
  {
    CGAL::Comparison_result operator()(Point_3, Point_3, Point_3)
    {return CGAL::ZERO;}
  };

  struct Has_on_bounded_side_3
  {
    //documented as Comparision_result
    CGAL::Comparison_result operator()(const Point_3&, const Point_3&, const Point_3&)
    {return CGAL::ZERO;}
    bool operator()(const Sphere_3&, const Point_3&)
    {return false;}
  };

  struct Compute_squared_radius_3
  {
    FT operator()(const Sphere_3&){return FT();}
  };

  struct Compute_squared_distance_3
  {
    FT operator()(const Point_3& , const Point_3& ){return FT();}
  };

  Compare_distance_3 compare_distance_3_object(){return Compare_distance_3();}
  Construct_sphere_3 construct_sphere_3_object(){return Construct_sphere_3();}
  Construct_projected_point_3 construct_projected_point_3_object(){return Construct_projected_point_3();}
  Compute_squared_distance_3 compute_squared_distance_3_object(){return Compute_squared_distance_3();}
  Do_intersect_3 do_intersect_3_object(){return Do_intersect_3();}
  Equal_3 equal_3_object(){return Equal_3();}
// } end of requirments from AABBGeomTraits


// requirements from SearchGeomTraits_3 {
  struct Iso_cuboid_3{};

  struct Construct_iso_cuboid_3{};

  struct Construct_min_vertex_3{};

  struct Construct_max_vertex_3{};

  struct Construct_center_3{};


  typedef const FT* Cartesian_const_iterator_3;


  struct Construct_cartesian_const_iterator_3{
    Construct_cartesian_const_iterator_3(){}
    Construct_cartesian_const_iterator_3(const Point_3&){}
    const FT* operator()(const Point_3&) const
    { return nullptr; }

    const FT* operator()(const Point_3&, int)  const
    { return nullptr; }
    typedef const FT* result_type;
  };
// } end of requirements from SearchGeomTraits_3

// requirements from SpatialSortingTraits_3 {
  struct Less_x_3{
    bool operator()(Point_3, Point_3){return false;}
  };
  struct Less_y_3{
    bool operator()(Point_3, Point_3){return false;}
  };
  struct Less_z_3{
    bool operator()(Point_3, Point_3){return false;}
  };
  Less_x_3 less_x_3_object()const{return Less_x_3();}
  Less_y_3 less_y_3_object()const{return Less_y_3();}
  Less_z_3 less_z_3_object()const{return Less_z_3();}
// } end of requirements from SpatialSortingTraits_3
};

namespace CGAL{

CGAL_DEFINE_COERCION_TRAITS_FOR_SELF(Custom_traits_Hausdorff::FT)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short, Custom_traits_Hausdorff::FT)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int, Custom_traits_Hausdorff::FT)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long, Custom_traits_Hausdorff::FT)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float, Custom_traits_Hausdorff::FT)
CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double, Custom_traits_Hausdorff::FT)

template<>struct Kernel_traits<Custom_traits_Hausdorff::Point_3>
{
  typedef Custom_traits_Hausdorff Kernel;
};

template<>
struct Algebraic_structure_traits< Custom_traits_Hausdorff::FT >
: public  Algebraic_structure_traits_base< Custom_traits_Hausdorff::FT, Field_tag >
{
};


template<>
struct Real_embeddable_traits< Custom_traits_Hausdorff::FT >
  : public INTERN_RET::Real_embeddable_traits_base< Custom_traits_Hausdorff::FT , CGAL::Tag_true>
{
  class To_double
    : public CGAL::cpp98::unary_function< Custom_traits_Hausdorff::FT, double > {
    public:
      double operator()( const Custom_traits_Hausdorff::FT&  ) const { return 0; }
  };
};

}//end CGAL


void exact(Custom_traits_Hausdorff::FT){}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

#include <CGAL/Polygon_mesh_processing/distance.h>

namespace PMP = CGAL::Polygon_mesh_processing;


typedef CGAL::Surface_mesh<K::Point_3> Mesh;

template <class GeomTraits, class TriangleMesh>
void general_tests(const TriangleMesh& m1,
                   const TriangleMesh& m2 )
{

  if (m1.number_of_vertices()==0 || m2.number_of_vertices()==0) return;

  std::vector<typename GeomTraits::Point_3> points;
    points.push_back(typename GeomTraits::Point_3(0,0,0));
    points.push_back(typename GeomTraits::Point_3(0,1,0));
    points.push_back(typename GeomTraits::Point_3(1,0,0));
    points.push_back(typename GeomTraits::Point_3(0,0,1));
    points.push_back(typename GeomTraits::Point_3(0,2,0));

  std::cout << "Symmetric distance between meshes (sequential) "
              << PMP::approximate_symmetric_Hausdorff_distance<CGAL::Sequential_tag>(
                  m1,m2,
                  PMP::parameters::number_of_points_per_area_unit(4000),
                  PMP::parameters::number_of_points_per_area_unit(4000))
              << "\n";

  std::cout << "Max distance to point set "
            << PMP::approximate_max_distance_to_point_set(m1,points,4000)
            << "\n";

  std::cout << "Max distance to triangle mesh (sequential) "
            << PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points,m1)
            << "\n";

  std::vector<typename GeomTraits::Point_3> samples;
  PMP::sample_triangle_mesh(m1, std::back_inserter(samples));
  std::cout << samples.size()<<" points sampled on mesh."<<std::endl;

}

void test_concept()
{
  typedef Custom_traits_Hausdorff CK;
  CGAL::Surface_mesh<CK::Point_3> m1, m2;
  general_tests<CK>(m1,m2);
}

int main(int argc, char** argv)
{
  if(argc != 3)
  {
    std::cerr << "Missing input meshes" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh m1,m2;
  std::ifstream input(argv[1]);
  input >> m1;
  input.close();
  input.open(argv[2]);
  input >> m2;
  input.close();
  std::cout << "First mesh has " << num_faces(m1) << " faces\n";
  std::cout << "Second mesh has " << num_faces(m2) << " faces\n";

  CGAL::Real_timer time;
  #if defined(CGAL_LINKED_WITH_TBB)
  time.start();
   std::cout << "Distance between meshes (parallel) "
             << PMP::approximate_Hausdorff_distance<CGAL::Parallel_tag>(
                  m1,m2,PMP::parameters::number_of_points_per_area_unit(4000))
             << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";
  #endif

  time.reset();
  time.start();
  std::cout << "Distance between meshes (sequential) "
            << PMP::approximate_Hausdorff_distance<CGAL::Sequential_tag>(
                 m1,m2,PMP::parameters::number_of_points_per_area_unit(4000))
            << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";

  general_tests<K>(m1,m2);

  test_concept();

  std::vector<std::vector<std::size_t> > faces;
  std::vector<K::Point_3> points;
  input.open(argv[1]);
  CGAL::IO::read_OFF(input, points, faces);
  input.close();

  std::vector<K::Point_3> samples;
  PMP::sample_triangle_soup(points, faces, std::back_inserter(samples));
  std::cout<<samples.size()<<" points sampled on soup."<<std::endl;


  return 0;
}
