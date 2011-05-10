
#include <CGAL/Pure_complex_data_structure.h>
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
    typedef typename TDS::Simplex Simplex;
    typedef typename TDS::Simplex_handle Simplex_handle;
    typedef typename TDS::Face Face;
    typedef typename TDS::Facet Facet;
    typedef typename TDS::Facet_iterator Facet_iterator;

    TDS pc(d);
    cerr << "\nChecking Pure_cds of (" << type << ") dimension "
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
        assert(!pc.is_simplex(Simplex_handle()));
        assert(pc.is_vertex(vhs[i]));
        assert(pc.is_simplex(vhs[i]->simplex()));

        if( pc.current_dimension() > 0 )
        {
            //int nbs = pc.number_of_simplices();
            pc.insert_in_simplex(pc.simplex(vhs[i+1]));
            ++nb_verts;
            //assert((size_t)(nbs+pc.current_dimension())==pc.number_of_simplices());
        }
        assert( pc.is_valid() );
    }
    assert((nb_verts == pc.number_of_vertices()));

	if( d > 1 )
    {
        // insert in hole
        std::vector<Simplex_handle> simps;
        simps.push_back(pc.simplices_begin());
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
            pc.gather_incident_upper_faces(vit, d, out);
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
        cerr << '\n' << pc.number_of_simplices() << " simplices, ";
        cerr << ' ' << nbfft << " facets.";
    }
 
    // TEST File I/O
    std::ofstream fo((string("output-pcds-")+type).c_str());
    if( d % 2 )
        CGAL::set_binary_mode(fo);
    fo << pc;
    fo.close();

    std::ifstream fi((string("output-pcds-")+type).c_str());
    if( d % 2 )
        CGAL::set_binary_mode(fi);
    TDS input_pcds(d);
    fi >> input_pcds;
    fi.close();

    // CLEAR
    pc.clear();
    assert(-2==pc.current_dimension());
    assert(pc.empty());
    assert( pc.is_valid() );
}

#define test_static(DIM) {  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dimension_tag<DIM> > TDS; \
        test<TDS>(DIM, string("static")+string(#DIM)); }
#define test_dyn(DIM) { \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dynamic_dimension_tag> TDS; \
        test<TDS>(DIM, string("dynamic")+string(#DIM)) ;}

#define test_mirror_static(DIM) {  \
    typedef CGAL::Pure_complex_ds_simplex<void, CGAL::TDS_simplex_mirror_storage_policy> My_ds_simplex;  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dimension_tag<DIM>, \
                                              CGAL::Pure_complex_ds_vertex<>, \
                                              My_ds_simplex> My_pcds;  \
        test<My_pcds>(DIM, string("mirror&static")+string(#DIM)); }
#define test_mirror_dyn(DIM) { \
    typedef CGAL::Pure_complex_ds_simplex<void, CGAL::TDS_simplex_mirror_storage_policy> My_ds_simplex;  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dynamic_dimension_tag, \
                                              CGAL::Pure_complex_ds_vertex<>, \
                                              My_ds_simplex> My_pcds;  \
        test<My_pcds>(DIM, string("mirror&dynamic")+string(#DIM)) ;}


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
