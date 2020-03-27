// Print a Triangulation<Traits> with 3d points in object file format (OFF).

#ifndef CGAL_TRIANGULATION_PRINT_OFF_H
#define CGAL_TRIANGULATION_PRINT_OFF_H 1

#include <map>

namespace CGAL {

template < class Triang >
void
triangulation_print_OFF( std::ostream& out,
                         const Triang& triang,
                         bool binary  = false,
                         bool noc     = false,
                         bool verbose = false) {
    CGAL_precondition( triang.is_valid());
    typedef typename Triang::Vertex           Vertex;
    typedef typename Triang::Vertex_iterator  Vertex_iterator;
    typedef typename Triang::Face_iterator    Face_iterator;
    // Build a map from vertex pointers to vertex indices.
    std::map<const Vertex*,std::size_t, std::less<const Vertex*> > mapping;
    std::size_t vn = 0;
    Vertex_iterator vi = triang.vertices_begin();
    for ( ; vi != triang.vertices_end(); ++vi) {
        CGAL_assertion( ! triang.is_infinite( vi));
        mapping[ &*vi] = vn;
        vn++;
    }
    CGAL_assertion( vn == std::size_t( triang.number_of_vertices()));

    // Count finite and infinite faces.
    std::size_t fn  = 0;
    Face_iterator fi = triang.faces_begin();
    for ( ; fi != triang.faces_end(); ++fi) {
        CGAL_assertion( ! triang.is_infinite( fi));
        fn++;
    }
    std::size_t fin = triang.number_of_faces() - fn;

    File_header_OFF  header( binary, noc, false, verbose);
    File_writer_OFF  writer( header);
    writer.write_header( out, vn, 3 * fn + fin, fn);

    vi = triang.vertices_begin();
    for ( ; vi != triang.vertices_end(); ++vi) {
        CGAL_assertion( ! triang.is_infinite( vi));
        writer.write_vertex( to_double(vi->point().x()),
                             to_double(vi->point().y()),
                             to_double(vi->point().z()));
    }
    writer.write_facet_header();

    fi = triang.faces_begin();
    while ( fn --) {
        writer.write_facet_begin( 3);
        CGAL_assertion( mapping.find(&*(fi->vertex(0))) != mapping.end());
        CGAL_assertion( mapping.find(&*(fi->vertex(1))) != mapping.end());
        CGAL_assertion( mapping.find(&*(fi->vertex(2))) != mapping.end());
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(0))]);
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(1))]);
        writer.write_facet_vertex_index( mapping[ &*(fi->vertex(2))]);
        writer.write_facet_end();
        ++fi;
    }
    CGAL_assertion( fi == triang.faces_end());
    writer.write_footer();
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_PRINT_OFF_H //
// EOF //
