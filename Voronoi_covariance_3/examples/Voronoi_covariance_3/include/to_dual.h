#ifndef TO_DUAL_H
#define TO_DUAL_H

#include <fstream>

// Convert a dual point to a point
template <typename R>
typename R::Point_3 to_dual (typename R::Plane_3 const& p) {
    typename R::Point_3 pp(-p.a() / p.d(), -p.b() / p.d(), -p.c() / p.d());

    return pp;
}

// Convert a dual plane to a plane
template <typename R>
CGAL::Point_triple<R> to_dual_plane (CGAL::Convex_hull_3::Plane_dual<R> const& p) {
    typename R::Plane_3 p1 = p.p1;
    typename R::Plane_3 p2 = p.p2;
    typename R::Plane_3 p3 = p.p3;

    CGAL::Point_triple<R> pp(to_dual<R>(p1),
                             to_dual<R>(p2),
                             to_dual<R>(p3));

    return pp;
}

// Write a dual polyhedron into an OFF file
template <typename R>
void convert_dual_OFF (std::string const& filename, Polyhedron_dual_3& P) {
    typedef Polyhedron_dual_3::Point_iterator                    Point_iterator;
    typedef Polyhedron_dual_3::Facet_iterator                    Facet_iterator;
    typedef Polyhedron_dual_3::Halfedge_around_facet_circulator  Halfedge_facet_circulator;

    std::ofstream file(filename.c_str());

    // 0. number of vertices / number of facets / number of edges
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
        << P.size_of_facets() << " 0" << std::endl;

    // 1. vertices definition
    for (Point_iterator vit = P.points_begin(); vit != P.points_end(); vit++) {
        typename R::Plane_3 p = *vit;
        file << to_dual<R>(p) << std::endl;
    }

    // 2. facets definition
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        file << std::endl;
    }

    file.close();
}

#endif

