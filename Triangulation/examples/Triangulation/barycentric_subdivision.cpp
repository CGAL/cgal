#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/internal/Combination_enumerator.h>
#include <iostream>
#include <vector>

using namespace std;

template< typename TDS >
void
make_face_from_vertices(
    const TDS & tds,
    const vector<typename TDS::Vertex_handle> & face_vertices,
    typename TDS::Face & face)
{ /* The main goal of this function is to find a full cell that contains a
     given set of vertices |face_vertices|. Then, it builds a corresponding
     |face|.
     */
    typedef typename TDS::Face                              Face;
    typedef typename TDS::Vertex_handle                     Vertex_handle;
    typedef typename TDS::Full_cell_handle                  Full_cell_handle;
    typedef typename TDS::Full_cell::Vertex_handle_iterator Vertex_h_iterator;
    
    // get the dimension of the face we want to build
    int fdim(face_vertices.size() - 1);
    if( fdim <= 0) exit(-1);
    
    // find all full cells incident to the first vertex of |face|
    typedef vector<Full_cell_handle> Cells;
    Cells cells;
    back_insert_iterator<Cells> out(cells);
    tds.incident_full_cells(face_vertices[0], out);
    // Iterate over the cells to find one which contains all the needed vertices
    for( typename Cells::iterator cit = cells.begin(); cit != cells.end(); ++cit )
    {
        // find if the cell *cit contains the Face |face|
        int i(0);
        for( ; i <= fdim; ++i )
        {
            Vertex_handle face_v(face_vertices[i]);
            bool found(false);
            for( Vertex_h_iterator vit = (*cit)->vertices_begin(); vit != (*cit)->vertices_end(); ++vit )
            {
                if( *vit == face_v )
                {
                    found = true;
                    break;
                }
            }
            if( ! found )
                break;
        }
        if( i > fdim ) // |*cit| contains |face|
        {
            face.set_full_cell(*cit);
            for( int i = 0; i <= fdim; ++i )
                face.set_index(i, (*cit)->index(face_vertices[i]));
            return;
        }
    }
    cerr << "\nCould not build a face from vertices\n";
    assert(false);
}

template< typename TDS >
void barycentric_subdivide(TDS & tds, typename TDS::Full_cell_handle fc)
{ /* This function builds the barycentric subdivision of a single full cell
     |fc| from a triangulation data structure |tds|.
     */
    typedef typename TDS::Full_cell_handle Full_cell_handle;
    typedef typename TDS::Vertex_handle Vertex_handle;
    typedef typename TDS::Face Face;
    
    const int dim = tds.current_dimension();

    // First, we read handles to the cell's vertices
    vector<Vertex_handle> vertices;
    vector<Vertex_handle> face_vertices;
    for( int i = 0; i <= dim; ++i )
        vertices.push_back(fc->vertex(i));

    // Then, we subdivide the cell |fc| once by inserting one vertex inside
    tds.insert_in_full_cell(fc);
    // From now on, we can't use the variable |fc|, since it has been subdivided.

    // Then, we subdivide facets of |fc| in order of decreasing dimension
    for( int d = dim-1; d > 0; --d )
    {
        face_vertices.resize(d+1);
        // create an enumerator of all faces of dimension d
        CGAL::internal::Combination_enumerator combi(d+1, 0, dim);
        while( ! combi.end() )
        {
            for( int i = 0; i <= d; ++i )
                face_vertices[i] = vertices[combi[i]];
            // we need to find a Full_cell that contains |face| as a... face.
            Face face(dim);
            make_face_from_vertices(tds, face_vertices, face);
            tds.insert_in_face(face);
            ++combi;
        }
    }
}

int main()
{
    const int sdim = 5; // dimension of TDS with compile-time dimension

    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<sdim> >
        TDS;
    typedef TDS::Vertex_handle    Vertex_handle;

    TDS  tds(sdim);

    Vertex_handle one_vertex;
    one_vertex = tds.insert_increase_dimension();
    assert( -1 == tds.current_dimension() );

    for( int i = 1; i < sdim+2; ++i )
        tds.insert_increase_dimension(one_vertex);
    assert( sdim   == tds.current_dimension() );
    assert( 2+sdim == tds.number_of_vertices() );
    assert( 2+sdim == tds.number_of_full_cells() );

    barycentric_subdivide(tds, tds.full_cells_begin());

    cout << "\nTriangulation has " << tds.number_of_full_cells() << " full cells";
    assert( tds.is_valid() );
    cout << " and is valid!\n";
    return 0;
}
