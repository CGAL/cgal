#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Iterator_project.h>

#include <cassert>

void test_range_based_loops(){

    typedef CGAL::Simple_cartesian<double>                     Kernel;
    typedef CGAL::Point_3<Kernel>                       Point;
    typedef CGAL::Plane_3<Kernel>                       Plane;
    typedef CGAL::Polyhedron_traits_3<Kernel>           Traits;

    typedef CGAL::Polyhedron_3<Traits>                  Polyhedron;

    typedef Polyhedron::HDS                     HDS;

    typedef Polyhedron::Vertex                  Vertex;
    CGAL_USE_TYPE(Polyhedron::Halfedge);
    CGAL_USE_TYPE(Polyhedron::Facet);

    typedef Polyhedron::Vertex_iterator         Vertex_iterator;
    typedef Polyhedron::Facet_iterator          Facet_iterator;
    typedef Polyhedron::Vertex_handle           Vertex_handle;
    typedef Polyhedron::Facet_handle            Facet_handle;
    typedef Polyhedron::Halfedge_handle         Halfedge_handle;
    typedef Polyhedron::Halfedge_iterator       Halfedge_iterator;


        Polyhedron P;
        Halfedge_handle h = P.make_tetrahedron(
                               Point( 1.0, 0.0, 0.0), Point( 0.0, 1.0, 0.0),
                               Point( 0.0, 0.0, 1.0), Point( 0.0, 0.0, 0.0));
        assert( P.is_valid());
        assert( P.is_tetrahedron( h));
    {
        typedef Polyhedron::Point_iterator Point_iterator;
        Point_iterator begin( P.points_begin());
        Point_iterator end( P.points_end());
        Vertex_iterator vitr = P.vertices_begin();
        for(Vertex_handle const &vh : P.vertex_handles()){
            assert( vh->point() == *begin);
            ++vitr;
            ++begin;
        }
        assert( vitr == P.vertices_end());
    }
    {
        typedef Polyhedron::Plane_iterator      Plane_iterator;
        Plane_iterator begin( P.planes_begin());
        Plane_iterator end( P.planes_end());
        Facet_iterator fitr = P.facets_begin();
        for(Facet_handle const &fh : P.facet_handles()){
            assert( fh->plane() == *begin);
            ++fitr;
            ++begin;
        }
        assert( fitr == P.facets_end());
    }
    {
        Halfedge_iterator hhitr = P.halfedges_begin();
        for(Halfedge_handle const & hh: P.halfedge_handles()){
            ++hhitr;
         }
        assert(hhitr == P.halfedges_end());
    }

}

int main() {
    test_range_based_loops();

    return 0;
}
