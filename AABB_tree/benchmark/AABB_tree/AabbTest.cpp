#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Timer.h>

#include "RaysGenerate.h"


typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

int main(int argc, char* argv[])
{
    const char* filename = (argc > 1) ? argv[1] : "data/data.ply";
    std::ifstream input(filename);
  
    Mesh mesh;
    CGAL::read_ply(input, mesh);

    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    Point p(0.0, 0.0, 0.0); /*POINT FOR SHOOTING RAY QUERIES*/

    int numberOfRays = 50000; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays); 
    CGAL::Timer time;
    time.start();

    for (int i=0; i<numberOfRays;i++){
        Vector v(rg.rayDirections[i]._x, rg.rayDirections[i]._y, rg.rayDirections[i]._z);
        Ray ray(p, v);
        
        Ray_intersection intersection = tree.first_intersection(ray);
        // if(intersection){
        //    if(boost::get<Point>(&(intersection->first))){
        //         const Point* p =  boost::get<Point>(&(intersection->first) );
        //         std::cout <<"Point of intersection : "<<  *p << std::endl;
        //     }
        // }
    }
    
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;   

    return 0;
}