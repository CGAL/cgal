#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

template <class TriangleMesh>
int genus(const TriangleMesh& tmesh)
{
    CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

    int euler_q = (int)CGAL::num_vertices(tmesh)
                - (int)CGAL::num_edges(tmesh)
                + (int)CGAL::num_faces(tmesh);

    int genus = (2 - euler_q) / 2;

    return genus;
}

bool test_on_file(const std::string& path)
{
    std::cout << std::endl;
    std::cout << "Path: " << path << std::endl;

    Mesh m;
    Polyhedron p;

    std::ifstream input(path);
    if (!input.is_open())
    {
        std::cerr << "Couldn't open the file" << std::endl;
        return false;
    }

    input >> m;

    input.close();
    input.clear();
    input.open(path);

    input >> p;

    int genus_m = genus(m);
    std::cout << "Genus of the mesh: " << genus_m << std::endl;

    int genus_p = genus(p);
    std::cout << "Genus of the polyhedron: " << genus_p << std::endl;

    assert(genus_m == genus_p);

    return true;
}

int main()
{
    // Test on tetrahedron
    {
        Polyhedron p;
        p.make_tetrahedron();

        int genus_p = genus(p);
        std::cout << "Genus of a tetrahedron: " << genus_p << std::endl;
        assert(genus_p == 0);
    }

    // Tests on various input files
    assert(test_on_file("data/elephant.off"));
    assert(test_on_file("data/cube.off"));
    assert(test_on_file("data/cheese.off"));
    assert(test_on_file("data/lion-head.off"));

    return EXIT_SUCCESS;
}

