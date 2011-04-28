#include <CGAL/Pure_complex_data_structure.h>
#include <iostream>
#include <vector>

int main()
{
    const int sdim = 27; // the dimension of the PCDS with compile-time dimension
    const int ddim = 15; // the dimension of the PCDS with dynamic dimension

    typedef CGAL::Pure_complex_data_structure<CGAL::Dynamic_dimension_tag>
            Dynamic_pcds;
    typedef CGAL::Pure_complex_data_structure<CGAL::Dimension_tag<sdim> >
            Static_pcds;
    typedef Static_pcds::Face           Face;
    typedef Static_pcds::Facet          Facet;
    typedef Static_pcds::Vertex_handle  Vertex_handle;
    typedef Static_pcds::Simplex_handle Simplex_handle;

    Static_pcds  S(158);  // the argument is not taken into account.
    Dynamic_pcds D(ddim); // the argument is taken into account.

    assert( sdim == S.ambient_dimension() );
    assert( ddim == D.ambient_dimension() );

    assert( -2 == S.current_dimension() );
    assert( S.is_valid() );

    std::vector<Vertex_handle> V(10);
    V[0] = S.insert_increase_dimension();
    assert( -1 == S.current_dimension() );

    for( int i = 1; i <= 5; ++i )
        V[i] = S.insert_increase_dimension(V[rand() % i]);

    assert( 4 == S.current_dimension() );
    assert( 6 == S.number_of_vertices() );
    assert( 6 == S.number_of_simplices() );

    Simplex_handle s = V[5]->simplex();

    V[6] = S.insert_in_simplex(s);
    
    assert( 7 == S.number_of_vertices() );
    assert( 10 == S.number_of_simplices() );

    s = V[3]->simplex();
    // ft will designate the Facet opposite to vertex 2 in s:
    Facet ft(s, 2);

    V[7] = S.insert_in_facet(ft);
    assert( 8 == S.number_of_vertices() );
    assert( 16 == S.number_of_simplices() );

    s = V[3]->simplex();
    // face will designate the edge joining vertices 2 and 4 in simplex s:
    Face face(s);
    face.set_index(0, 2);
    face.set_index(1, 4);

    V[8] = S.insert_in_face(face);
    assert( S.is_valid() );

    Simplex_handle hole[2];
    hole[0] = V[8]->simplex();
    hole[1] = hole[0]->neighbor(0);
    // ft will be a face on the boundary of the union of hole[0] and hole[1]:
    ft = Facet(hole[0], 1);

    V[9] = S.insert_in_hole(hole, hole+2, ft);
    assert( S.is_valid() );

    std::cout << "Ok.\n";
    return 0;
}
