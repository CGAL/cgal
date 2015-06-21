#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/OFF_to_nef_3.h>

int main()
{
   typedef CGAL::Homogeneous<CGAL::Exact_integer>  Kernel;
   typedef CGAL::Nef_polyhedron_3<Kernel> Nef_3;

   Nef_3 N;
   std::size_t discarded = CGAL::OFF_to_nef_3 (std::cin, N, true);

   std::cout << "Nef vertices: "
             << N.number_of_vertices() << std::endl;
   std::cout << "Nef edges: "
             << N.number_of_edges() << std::endl;
   std::cout << "Nef facets: "
             << N.number_of_facets() << std::endl;
   std::cout << "Nef volumes: "
             << N.number_of_volumes() << std::endl;
   std::cout << "number of discarded facets: "
             << discarded << std::endl;

   return 0;
}
