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

  void progress_filtering_intersections(double d)
   {
    d /= normalize;
    total += d;
    if(total > bound){
      std::cout << std::setprecision(3) << total*100 << " %   in " << std::setprecision(5) << t.time() << " sec." << std::endl;
      bound += 0.1;
    }
  }

  void start_triangulating_faces(std::size_t tf)
  {
    tfaces = tf;
    bound_faces = tf/10;
  }

  void face_triangulation(std::size_t i)
  {
    if(i> bound_faces){
      std::cout << double(i)/double(tfaces) * 100  << " %" << std::endl;
      bound_faces += tfaces/10;
    }
  }

  void start_coplanar_faces(std::size_t tc)
  {
    std::cout << "Visitor::start_coplanar_faces() at " << t.time() << " sec." << std::endl;
    tcoplanar= tc;
    count_coplanar = 0;
    bound_coplanar = tcoplanar/10;
  }

  void intersection_of_coplanar_faces_step()
  {
    ++count_coplanar;
    if(count_coplanar> bound_coplanar){
      std::cout << "Visitor::coplanar_faces: " << double(count_coplanar)/double(tcoplanar) * 100  << " % " << std::endl;
      bound_coplanar += tcoplanar/10;
    }
  }

  void start_intersection_points(std::size_t ti)
  {
    std::cout << "Visitor::start_intersection_points() at " << t.time() << " sec." << std::endl;
    tintersection= ti;
    count_intersection = 0;
    bound_intersection = tintersection/10;
  }

  void edge_face_intersections_step()
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
  std::size_t count = 0;

  std::size_t bound_faces = 0;
  std::size_t tfaces = 0;

  std::size_t bound_coplanar = 0;
  std::size_t tcoplanar = 0;
  std::size_t count_coplanar = 0;

  std::size_t bound_intersection = 0;
  std::size_t tintersection = 0;
  std::size_t count_intersection = 0;
  CGAL::Timer t;
};


struct Visitor :
  public PMP::Corefinement::Default_visitor<Mesh>
{
  std::shared_ptr<Visitor_rep> sptr;
  mutable std::size_t tf_counter = 0;

  Visitor()
    : sptr(std::make_shared<Visitor_rep>())
  {}

  void progress_filtering_intersections(double d)
  {
    sptr->progress_filtering_intersections(d);
  }

  void start_filtering_intersections() const
  {
    std::cout << "Visitor::start_filtering_intersections() at " << sptr->time() << " sec." << std::endl;
  }
  void end_filtering_intersections() const
  {
    std::cout << "Visitor::end_filtering_intersections() at " << sptr->time() << " sec."  << std::endl;
  }

  void start_triangulating_faces(std::size_t tf) const
  {
    std::cout << "Visitor::start_triangulation() with " << tf << " faces at " << sptr->time() << " sec."  << std::endl;
    sptr->start_triangulating_faces(tf);
    tf_counter = 0;
  }

  void triangulating_faces_step() const
  {
    sptr->face_triangulation(tf_counter++);
  }

  void end_triangulating_faces()const
  {
    std::cout << "Visitor::end_triangulating_faces() at " << sptr->time() << " sec."  << std::endl;
  }

  void start_handling_intersection_of_coplanar_faces(std::size_t i) const
  {
    sptr->start_coplanar_faces(i);
  }

  void intersection_of_coplanar_faces_step() const
  {
    sptr->intersection_of_coplanar_faces_step();
  }

  void end_handling_intersection_of_coplanar_faces() const
  {
    std::cout << "Visitor::end_coplanar_faces() at " << sptr->time() << " sec." << std::endl;
  }

  void start_handling_edge_face_intersections(std::size_t i) const
  {
    sptr->start_intersection_points(i);
  }

  void edge_face_intersections_step() const
  {
    sptr->edge_face_intersections_step();
  }

  void end_handling_edge_face_intersections() const
  {
    std::cout << "Visitor::end_intersection_points() at " << sptr->time() << " sec." << std::endl;
  }

  void start_building_output() const
  {
    std::cout << "Visitor::start_building_output() at " << sptr->time()  << " sec."<< std::endl;
  }

  void end_building_output() const
  {
    std::cout << "Visitor::end_building_output() at " << sptr->time() << " sec." << std::endl;
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
