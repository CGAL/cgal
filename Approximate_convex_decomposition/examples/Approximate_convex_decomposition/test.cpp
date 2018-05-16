#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>
#include <fstream>
#include <string>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

template <typename size_type>
int genus_aux(size_type vertices_cnt, size_type edges_cnt, size_type faces_cnt)
{
    int euler_q = (int)vertices_cnt
                - (int)edges_cnt
                + (int)faces_cnt;

    int genus = (2 - euler_q) / 2;

    return genus;
}

int mesh_genus(const Mesh& mesh)
{
    return genus_aux<Mesh::size_type>(mesh.number_of_vertices(),
                                      mesh.number_of_edges(),
                                      mesh.number_of_faces());
}

int polyhedron_genus(const Polyhedron& polyhedron)
{
    return genus_aux<Polyhedron::size_type>(polyhedron.size_of_vertices(),
                                            polyhedron.size_of_halfedges() / 2,
                                            polyhedron.size_of_facets());
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

    if (!CGAL::is_triangle_mesh(p) ||
        !CGAL::is_triangle_mesh(m))
    {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return false;
    }

    int genus_m = mesh_genus(m);
    std::cout << "Genus of the mesh: " << genus_m << std::endl;

    int genus_p = polyhedron_genus(p);
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

        int genus = polyhedron_genus(p);
        std::cout << "Genus of a tetrahedron: " << genus << std::endl;
        assert(genus == 0);
    }

    // Tests on various input files
    assert(test_on_file("data/elephant.off"));
    assert(test_on_file("data/cube.off"));
    assert(test_on_file("data/cheese.off"));
    assert(test_on_file("data/lion-head.off"));

    return EXIT_SUCCESS;
}
