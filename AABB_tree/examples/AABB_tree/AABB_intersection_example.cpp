#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_trees/intersection.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;

typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive_3<K, Iterator> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);

    Point e(1.2, 0.2, 0.2);
    Point f(0.2, 1.2, 0.2);
    Point g(0.2, 0.2, 1.2);
    Point h(0.2, 0.2, 0.2);

    std::vector<Triangle> tetra1;
    tetra1.push_back(Triangle(a,b,c));
    tetra1.push_back(Triangle(a,b,d));
    tetra1.push_back(Triangle(a,d,c));
    tetra1.push_back(Triangle(b,c,d));

    std::vector<Triangle> tetra2;
    tetra2.push_back(Triangle(e,f,g));
    tetra2.push_back(Triangle(e,f,h));
    tetra2.push_back(Triangle(e,h,g));
    tetra2.push_back(Triangle(f,g,h));

    // constructs AABB tree
    Tree tree1(tetra1.begin(),tetra1.end());
    Tree tree2(tetra2.begin(),tetra2.end());

    // do intersect
    std::cout << "Tetrahedrons do intersect: " << CGAL::AABB_trees::do_intersect(tree1, tree2) << std::endl;

    std::vector< std::pair<Iterator, Iterator> > intersections;
    CGAL::AABB_trees::all_pairs_of_intersecting_primitives(tree1, tree2, std::back_inserter(intersections));

    for(auto [id1, id2]: intersections)
      std::cout << "ids: " << std::distance(tetra1.begin(), id1) << " " << std::distance(tetra2.begin(), id2) << std::endl;

    return EXIT_SUCCESS;
}
