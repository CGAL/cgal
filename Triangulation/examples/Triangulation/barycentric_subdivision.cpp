#include <CGAL/Triangulation_data_structure.h>
#include <CGAL/internal/Combination_enumerator.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>

template< typename TDS >
void find_face_from_vertices(const TDS & tds,
        const std::vector<typename TDS::Vertex_handle> & face_vertices,
        typename TDS::Face & face);

template< typename TDS >
void barycentric_subdivide(TDS & tds, typename TDS::Full_cell_handle fc)
{ /* This function builds the barycentric subdivision of a single
     full cell |fc| from a triangulation data structure |tds|.  */
    typedef typename TDS::Vertex_handle Vertex_handle;
    typedef typename TDS::Face Face;
    const int dim = tds.current_dimension();

    // First, read handles to the cell's vertices
    std::vector<Vertex_handle> vertices;
    std::vector<Vertex_handle> face_vertices;
    for( int i = 0; i <= dim; ++i ) vertices.push_back(fc->vertex(i));

    // Then, subdivide the cell |fc| once by inserting one vertex
    tds.insert_in_full_cell(fc);
    // From now on, we can't use the variable |fc|...

    // Then, subdivide faces of |fc| in order of decreasing dimension
    for( int d = dim-1; d > 0; --d )
    {
        face_vertices.resize(d+1);
        // The following class
        // enumerates all (d+1)-tuple of the set {0, 1, ..., dim}
        CGAL::internal::Combination_enumerator combi(d+1, 0, dim);
        while( ! combi.end() )
        {
            for( int i = 0; i <= d; ++i )
                face_vertices[i] = vertices[combi[i]];
            // we need to find a face with face_vertices
            Face face(dim);
            find_face_from_vertices(tds, face_vertices, face);
            tds.insert_in_face(face);
            ++combi;
        }
    }
}

template< typename TDS >
void find_face_from_vertices( const TDS & tds,
        const std::vector<typename TDS::Vertex_handle> & face_vertices,
        typename TDS::Face & face)
{ /* The main goal of this function is to find a full cell that
     contains a given set of vertices |face_vertices|. Then, it
     builds a corresponding |face|. */
    typedef typename TDS::Vertex_handle           Vertex_handle;
    typedef typename TDS::Full_cell_handle        Full_cell_handle;
    typedef typename TDS::Full_cell::Vertex_handle_iterator Vertex_h_iterator;

    // get the dimension of the face we want to build
    std::size_t fdim(face_vertices.size() - 1);
    if( fdim <= 0) exit(-1);

    // find all full cells incident to the first vertex of |face|
    typedef std::vector<Full_cell_handle> Cells;
    Cells cells;
    std::back_insert_iterator<Cells> out(cells);
    tds.incident_full_cells(face_vertices[0], out);
    // Iterate over the cells to find one which contains the face_vertices
    for( typename Cells::iterator cit = cells.begin(); cit != cells.end(); ++cit){
        // find if the cell *cit contains the Face |face|
        std::size_t i = 0;
        for( ; i <= fdim; ++i ) {
            Vertex_handle face_v(face_vertices[i]);
            bool found(false);
            Vertex_h_iterator vit = (*cit)->vertices_begin();
            for( ; vit != (*cit)->vertices_end(); ++vit ) {
                if( *vit == face_v ) {
                    found = true;
                    break;
                }
            }
            if( ! found )
                break;
        }
        if( i > fdim ) {// the full cell |*cit| contains |face|
            face.set_full_cell(*cit);
            for( std::size_t i = 0; i <= fdim; ++i )
            {
              face.set_index(static_cast<int>(i), 
                             (*cit)->index(face_vertices[i]));
            }
            return;
        }
    }
    std::cerr << "Could not build a face from vertices"<<std::endl;
    CGAL_assertion(false);
}


int main()
{
    const int sdim = 5; // dimension of TDS with compile-time dimension
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<sdim> >
        TDS;
    TDS  tds(sdim);

    TDS::Vertex_handle one_vertex = tds.insert_increase_dimension();
    for( int i = 1; i < sdim+2; ++i )
        tds.insert_increase_dimension(one_vertex);
    // we get a triangulation of space of dim sdim homeomorphic to
    // the boundary of simplex of dimension sdim+1 with sdim+2 vertices
    CGAL_assertion( sdim   == tds.current_dimension() );
    CGAL_assertion( 2+sdim == tds.number_of_vertices() );
    CGAL_assertion( 2+sdim == tds.number_of_full_cells() );

    barycentric_subdivide(tds, tds.full_cells_begin());

    // The number of full cells should be twice the factorial of
    // |tds.current_dimension()+1|. Eg, 1440 for dimension 5.
    std::cout << "Triangulation has " 
        << tds.number_of_full_cells() << " full cells";
    CGAL_assertion( tds.is_valid() );
    std::cout << " and is valid!"<<std::endl;
    return 0;
}
