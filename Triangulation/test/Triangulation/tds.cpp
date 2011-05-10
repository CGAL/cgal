
#include <CGAL/Triangulation_data_structure.h>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

template<typename TDS>
void test(const int d, const string & type)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that TDS contains sub-types, before
    // instanciating the function.
    typedef typename TDS::Vertex Vertex;
    typedef typename TDS::Vertex_handle Vertex_handle;
    typedef typename TDS::Vertex_iterator Vertex_iterator;
    typedef typename TDS::Full_cell Full_cell;
    typedef typename TDS::Full_cell_handle Full_cell_handle;
    typedef typename TDS::Face Face;
    typedef typename TDS::Facet Facet;
    typedef typename TDS::Facet_iterator Facet_iterator;

    TDS pc(d);
    cerr << "\nChecking Tds of (" << type << ") dimension "
        << pc.ambient_dimension();
    assert(pc.empty());
    vector<Vertex_handle> vhs;
    vhs.push_back(pc.insert_increase_dimension());
    assert( pc.is_valid() );
    size_t nb_verts = 1;
    for( int i = 0; i <= d; ++i )
    {
        vhs.push_back(pc.insert_increase_dimension(vhs[0]));
        ++nb_verts;
        
        assert(i == pc.current_dimension());
        assert(!pc.is_vertex(Vertex_handle()));
        assert(!pc.is_full_cell(Full_cell_handle()));
        assert(pc.is_vertex(vhs[i]));
        assert(pc.is_full_cell(vhs[i]->full_cell()));

        if( pc.current_dimension() > 0 )
        {
            //int nbs = pc.number_of_full_cells();
            pc.insert_in_full_cell(pc.full_cell(vhs[i+1]));
            ++nb_verts;
            //assert((size_t)(nbs+pc.current_dimension())==pc.number_of_full_cells());
        }
        assert( pc.is_valid() );
    }
    assert((nb_verts == pc.number_of_vertices()));

	if( d > 1 )
    {
        // insert in hole
        std::vector<Full_cell_handle> simps;
        simps.push_back(pc.full_cells_begin());
        simps.push_back(pc.neighbor(simps[0],0));
        pc.insert_in_hole(simps.begin(), simps.end(), Facet(simps[0],1));
    }

   // TEST Faces enumeration
    typedef std::vector<Face> Faces;
    Faces faces;
    for( Vertex_iterator vit = pc.vertices_begin(); vit != pc.vertices_end(); ++vit )
    {
        for( int d = 1; d < pc.current_dimension() - 1; ++d )
        {
            cerr << '\n' << d << "-dimensional faces adjacent to " << &(*vit)
                << " ( current dimension is " << pc.current_dimension() << " )";
            faces.clear();
            std::back_insert_iterator<Faces> out(faces);
            pc.incident_upper_faces(vit, d, out);
            typename Faces::iterator fit = faces.begin();
            while( fit != faces.end() )
            {
                cerr << '\n';
                for( int i = 0; i <= d; ++i )
                    cerr << ' ' << &(*fit->vertex(i));
                ++fit;
            }
        }
    }

    // TEST Finite iterators
    if( pc.current_dimension() > 0 )
    {
        Facet_iterator fit = pc.facets_begin();
        size_t nbfft(0);
        while( fit != pc.facets_end() )
            ++fit, ++nbfft;
        cerr << '\n' << pc.number_of_full_cells() << " full cells, ";
        cerr << ' ' << nbfft << " facets.";
    }
 
    // TEST File I/O
    std::ofstream fo((string("output-tds-")+type).c_str());
    if( d % 2 )
        CGAL::set_binary_mode(fo);
    fo << pc;
    fo.close();

    std::ifstream fi((string("output-tds-")+type).c_str());
    if( d % 2 )
        CGAL::set_binary_mode(fi);
    TDS input_tds(d);
    fi >> input_tds;
    fi.close();

    // CLEAR
    pc.clear();
    assert(-2==pc.current_dimension());
    assert(pc.empty());
    assert( pc.is_valid() );
}

#define test_static(DIM) {  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<DIM> > TDS; \
        test<TDS>(DIM, string("static")+string(#DIM)); }
#define test_dyn(DIM) { \
    typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag> TDS; \
        test<TDS>(DIM, string("dynamic")+string(#DIM)) ;}

#define test_mirror_static(DIM) {  \
    typedef CGAL::Triangulation_ds_full_cell<void, CGAL::TDS_full_cell_mirror_storage_policy> My_ds_full_cell;  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<DIM>, \
                                              CGAL::Triangulation_ds_vertex<>, \
                                              My_ds_full_cell> My_tds;  \
        test<My_tds>(DIM, string("mirror&static")+string(#DIM)); }
#define test_mirror_dyn(DIM) { \
    typedef CGAL::Triangulation_ds_full_cell<void, CGAL::TDS_full_cell_mirror_storage_policy> My_ds_full_cell;  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag, \
                                              CGAL::Triangulation_ds_vertex<>, \
                                              My_ds_full_cell> My_tds;  \
        test<My_tds>(DIM, string("mirror&dynamic")+string(#DIM)) ;}


int main()
{
    //test_static(11);
    //test_static(10);
    test_static(4);
    test_static(3);
    test_static(2);
    test_static(1);

    //test_dyn(11);
    //test_dyn(10);
    test_dyn(4);
    test_dyn(3);
    test_dyn(2);
    test_dyn(1);

    //test_mirror_static(11);
    //test_mirror_static(10);
    test_mirror_static(4);
    test_mirror_static(3);
    test_mirror_static(2);
    test_mirror_static(1);

    //test_mirror_dyn(11);
    //test_mirror_dyn(10);
    test_mirror_dyn(4);
    test_mirror_dyn(3);
    test_mirror_dyn(2);
    test_mirror_dyn(1);

    cerr << std::endl;
    return 0;
}
