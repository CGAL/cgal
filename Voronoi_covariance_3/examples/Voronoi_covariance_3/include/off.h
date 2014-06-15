#ifndef _OFF_H_
#define _OFF_H_

#include <fstream>

// Write into an OFF file to visualize with GeomView
template <typename K, typename Polyhedron>
void convertToOFF (std::string const& filename, Polyhedron& P) {
    typedef typename Polyhedron::Facet_iterator                   Facet_iterator;
    typedef typename Polyhedron::Halfedge_handle                   Halfedge_handle;
    typedef typename Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

    std::ofstream file(filename.c_str());

    // 0. number of vertices / number of facets / number of edges
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
        << P.size_of_facets() << " 0" << std::endl;

    // 1. vertices definition
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<typename K::Point_3>( file, "\n"));

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

