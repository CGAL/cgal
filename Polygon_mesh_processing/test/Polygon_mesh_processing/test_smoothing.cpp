#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/property_map/property_map.hpp>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

void calc_angles(Mesh& mesh, double& min_a, double& max_a, double& avg_a)
{
    using namespace boost::accumulators;
    accumulator_set< double,
      features< tag::min, tag::max, tag::mean > > acc;

    typename boost::property_map<Mesh, CGAL::vertex_point_t>::type
      vpmap = get(CGAL::vertex_point, mesh);

    typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
    double rad_to_deg = 180. / CGAL_PI;

    BOOST_FOREACH(halfedge_descriptor h, halfedges(mesh))
    {
        if(face(h, mesh) == boost::graph_traits<Mesh>::null_face())
            continue;

        typename K::Point_3 a = get(vpmap, source(h, mesh));
        typename K::Point_3 b = get(vpmap, target(h, mesh));
        typename K::Point_3 c = get(vpmap, target(next(h, mesh), mesh));
        typename K::Vector_3 ba(b, a);
        typename K::Vector_3 bc(b, c);

        double cos_angle = (ba * bc)
          / std::sqrt(ba.squared_length() * bc.squared_length());

        acc(std::acos(cos_angle) * rad_to_deg);
    }

    min_a = extract_result< tag::min >(acc);
    max_a = extract_result< tag::max >(acc);
    avg_a = extract_result< tag::mean >(acc);
}

void calc_areas(Mesh& mesh, double& min_a, double& max_a, double& avg_a)
{
    using namespace boost::accumulators;
    accumulator_set< double,
      features< tag::min, tag::max, tag::mean > > acc;

    typename boost::property_map<Mesh, CGAL::vertex_point_t>::type
      vpmap = get(CGAL::vertex_point, mesh);

    typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;

    BOOST_FOREACH(face_descriptor f, faces(mesh))
    {
        if(f == boost::graph_traits<Mesh>::null_face())
            continue;

        double area = CGAL::Polygon_mesh_processing::face_area(f, mesh);
        acc(area);
    }

    min_a = extract_result< tag::min >(acc);
    max_a = extract_result< tag::max >(acc);
    avg_a = extract_result< tag::mean >(acc);
}

template<class T>
bool check_value_equal(T a, T b)
{
  if (abs(a-b) > 1e-3)
  {
    std::cerr << "Value not equal! ";
    std::cerr << a << " is not " << b << std::endl;
    return false;
  }
  return true;
}


int main(int argc, char* argv[]){

    std::string filename;
    std::ifstream input;
    Mesh mesh;
    double min_a, max_a, mean_a;

    filename = "data/curved_polygon.off";
    input.open(filename);
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid .off file." << std::endl;
        return 1;
    }
    input.close();

    calc_angles(mesh, min_a, max_a, mean_a);

    CGAL::Polygon_mesh_processing::smooth_angles(mesh);
    calc_angles(mesh, min_a, max_a, mean_a);

    if(!check_value_equal(min_a, 24.980))
    {
        return EXIT_FAILURE;
    }
    if(!check_value_equal(max_a, 108.723))
    {
        return EXIT_FAILURE;
    }

    input.open(filename);
    input>>mesh;
    input.close();

    CGAL::Polygon_mesh_processing::smooth_areas(mesh);
    calc_areas(mesh, min_a, max_a, mean_a);

    if(!check_value_equal(min_a, 0.476))
    {
        return EXIT_FAILURE;
    }
    if(!check_value_equal(max_a, 0.567))
    {
        return EXIT_FAILURE;
    }
    if(!check_value_equal(mean_a, 0.542))
    {
        return EXIT_FAILURE;
    }

    return 0;
}
