#include <iostream>
#include <fstream>
#include <cstdlib>

#include <CGAL/Surface_mesh_segmentation.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Surface_mesh_segmentation<Polyhedron> Segmentation;

int main(int argc, char **argv)
{
    std::ifstream meshFile("camel.off");
    if(!meshFile) 
    { 
        std::cerr << "Could not open the file!" << std::endl;
        return EXIT_FAILURE;
    }
    
    Polyhedron mesh;   
    meshFile >> mesh;
    
    int ray_sqrt = 7; // cast (7*7)49 rays per facet
    int number_of_clusters = 2; // apply gmm fitting with 2 clusters
    double cone_angle =  (2.0 / 3.0) * CGAL_PI;

    Segmentation segmentation(&mesh, ray_sqrt, cone_angle, number_of_clusters);
    //segmentation.write_sdf_values("sdf_values.txt");
}