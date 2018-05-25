#include <CGAL/approx_decomposition.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

int main()
{
    // reading mesh
    Polyhedron mesh;
    
    std::ifstream input("data/elephant.off");
    
    if (!input || !(input >> mesh))
    {
        std::cout << "Failed to read mesh" << std::endl;
        return EXIT_FAILURE;
    }

    if (mesh.empty() || !CGAL::is_triangle_mesh(mesh) || !CGAL::is_closed(mesh))
    {
        std::cout << "Input mesh is invalid" << std::endl;
        return EXIT_FAILURE;
    }

    // computing concavity value
    //double concavity = CGAL::concavity_value(mesh);

    // writing result
    //std::cout << concavity << std::endl;

    return EXIT_SUCCESS;
}

