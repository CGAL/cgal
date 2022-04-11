#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Visitor_rep{

  Visitor_rep(double normalize = 4)
    : normalize(normalize)
  {
    t.start();
  }

  void progress_filter_intersection(double d)
   {
    d /= normalize;
    total += d;
    if(total > bound){
      std::cout << std::setprecision(3) << total*100 << " %   in " << std::setprecision(5) << t.time() << " sec." << std::endl;
      bound += 0.1;
    }
  }

  void start_triangulation(int tf)
  {
    tfaces = tf;
    bound_faces = tf/10;
  }

  void progress_triangulation(int i)
  {
    if(i> bound_faces){
      std::cout << double(i)/double(tfaces) * 100  << " %" << std::endl;
      bound_faces += tfaces/10;
    }
  }

  void start_coplanar_faces(int tc)
  {
    std::cout << "Visitor::start_coplanar_faces() at " << t.time() << " sec." << std::endl;
    tcoplanar= tc;
    count_coplanar = 0;
    bound_coplanar = tcoplanar/10;
  }

  void coplanar_faces_step()
  {
    ++count_coplanar;
    if(count_coplanar> bound_coplanar){
      std::cout << "Visitor::coplanar_faces: " << double(count_coplanar)/double(tcoplanar) * 100  << " % " << std::endl;
      bound_coplanar += tcoplanar/10;
    }
  }

  void start_intersection_points(int ti)
  {
    std::cout << "Visitor::start_intersection_points() at " << t.time() << " sec." << std::endl;
    tintersection= ti;
    count_intersection = 0;
    bound_intersection = tintersection/10;
  }

  void intersection_points_step()
  {
    ++count_intersection;
    if(count_intersection> bound_intersection){
      std::cout << "Visitor::intersection_points: " << double(count_intersection)/double(tintersection) * 100  << " % " << std::endl;
      bound_intersection += tintersection/10;
    }
  }

  double time() const
  {
    return t.time();
  }

  double normalize;
  double bound = 0.1;
  double total = 0;
  int count = 0;

  int bound_faces = 0;
  int tfaces = 0;

  int bound_coplanar = 0;
  int tcoplanar = 0;
  int count_coplanar = 0;

  int bound_intersection = 0;
  int tintersection = 0;
  int count_intersection = 0;
  CGAL::Timer t;
};


struct Visitor :
  public PMP::Corefinement::Default_visitor<Mesh>
{
  std::shared_ptr<Visitor_rep> sptr;

  Visitor()
    : sptr(std::make_shared<Visitor_rep>())
  {}

  void progress_filter_intersection(double d)
  {
    sptr->progress_filter_intersection(d);
  }

  void start_filter_intersections() const
  {
    std::cout << "Visitor::start_filter_intersections() at " << sptr->time() << " sec." << std::endl;
  }
  void end_filter_intersections() const
  {
    std::cout << "Visitor::end_filter_intersections() at " << sptr->time() << " sec."  << std::endl;
  }

  void start_triangulation(int tf)const
  {
    std::cout << "Visitor::start_triangulation() with " << tf << " faces at " << sptr->time() << " sec."  << std::endl;
    sptr->start_triangulation(tf);
  }

  void progress_triangulation(int i) const
  {
    sptr->progress_triangulation(i);
  }

  void end_triangulation()const
  {
    std::cout << "Visitor::end_triangulation() at " << sptr->time() << " sec."  << std::endl;
  }

  void start_coplanar_faces(int i) const
  {
    sptr->start_coplanar_faces(i);
  }

  void coplanar_faces_step() const
  {
    sptr->coplanar_faces_step();
  }

  void end_coplanar_faces() const
  {
    std::cout << "Visitor::end_coplanar_faces() at " << sptr->time() << " sec." << std::endl;
  }

  void start_intersection_points(int i) const
  {
    sptr->start_intersection_points(i);
  }

  void intersection_points_step() const
  {
    sptr->intersection_points_step();
  }

  void end_intersection_points() const
  {
    std::cout << "Visitor::end_intersection_points() at " << sptr->time() << " sec." << std::endl;
  }

  void start_build_output() const
  {
    std::cout << "Visitor::start_build_output() at " << sptr->time()  << " sec."<< std::endl;
  }

  void build_output_step() const
  {
    std::cout << "Visitor::build_output_step() at " << sptr->time() << " sec." << std::endl;
  }

  void end_build_output() const
  {
    std::cout << "Visitor::end_build_output() at " << sptr->time() << " sec." << std::endl;
  }
};


int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh1, mesh2;
  if(!CGAL::IO::read_polygon_mesh(filename1, mesh1) || !CGAL::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Timer t;
  t.start();
  Mesh out;
  Visitor visitor;

  bool valid_union = PMP::corefine_and_compute_union (mesh1,mesh2, out, CGAL::parameters::visitor(visitor));

  std::cout << "Global timer = " << t.time() << " sec." << std::endl;


  if(valid_union)
  {
    std::cout << "Union was successfully computed\n";
    CGAL::IO::write_polygon_mesh("union.off", out, CGAL::parameters::stream_precision(17));
    return 0;
  }

  std::cout << "Union could not be computed\n";

  return 1;
}
