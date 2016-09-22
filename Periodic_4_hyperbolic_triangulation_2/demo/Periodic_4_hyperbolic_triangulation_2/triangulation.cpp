

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

#include <CGAL/Bit_word_utility_4.h>

typedef CORE::Expr                                              					NT;					
typedef CGAL::Cartesian<NT>                                     					Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> 		Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>				Triangulation;

typedef Triangulation::Point_2	                                   					Point;
typedef Triangulation::Face_handle													Face_handle;
typedef Triangulation::Vertex_handle												Vertex_handle;
typedef Triangulation::Locate_type 													Locate_type;
typedef Triangulation::Edge                                                         Edge;
typedef Traits::FT                                                                  FT;
typedef Triangulation::Face_iterator                                                Face_iterator;
typedef Triangulation::Offset                                                       Offset;
typedef Traits::Side_of_fundamental_octagon                                         Side_of_fundamental_octagon;

typedef unsigned short int                                                          Int;



using namespace CGAL;

template<class GT>
Hyperbolic_octagon_translation_matrix<GT> get_matrix(vector<Int> seq) {
    Hyperbolic_octagon_translation_matrix<GT> r;
    vector<Hyperbolic_octagon_translation_matrix<GT> > gens;
    get_generators(gens);
    for (int i = 0; i < seq.size(); i++) {
        r = r * gens[seq[i]];
    }
    return r;
}

void get_offsets(vector<Offset>& v) {
    for (int i = 0; i < 8; i++) {
        v.push_back(Offset(i));
    }
    for (int i = 1; i < 4; i++) { // word length
        vector<Offset> old;
        for (int j = 0; j < v.size(); j++) {
            old.push_back(v[j]);
        }
        for (int j = 0; j < old.size(); j++) {
            for (int k = 0; k < 8; k++) {
                v.push_back(old[j].append(Offset(k)));
            }
        }
    }
}


int main(void) {

    //Offset o, s;
    //vector<Int> v, r;

    vector<Offset> ov;
    
    ov.push_back( Offset(4, 1, 6, 3) );
    ov.push_back( Offset(4, 1, 6) );
    ov.push_back( Offset(4, 1) );
    ov.push_back( Offset(4) );
    ov.push_back( Offset(4, 7) );
    ov.push_back( Offset(4, 7, 2) );
    ov.push_back( Offset(5, 2, 7, 4) );
    ov.push_back( Offset(5, 2, 7) );
    ov.push_back( Offset(5, 2) );
    ov.push_back( Offset(5) );
    ov.push_back( Offset(5, 0) );
    ov.push_back( Offset(5, 0, 3) );
    ov.push_back( Offset(6, 3, 0, 5) );
    ov.push_back( Offset(6, 3, 0) );
    ov.push_back( Offset(6, 3) );
    ov.push_back( Offset(6) );
    ov.push_back( Offset(6, 1) );
    ov.push_back( Offset(6, 1, 4) );
    ov.push_back( Offset(7, 4, 1, 6) );
    ov.push_back( Offset(7, 4, 1) );
    ov.push_back( Offset(7, 4) );
    ov.push_back( Offset(7) );
    ov.push_back( Offset(7, 2) );
    ov.push_back( Offset(7, 2, 5) );
    ov.push_back( Offset(0, 5, 2, 7) );
    ov.push_back( Offset(0, 5, 2) );
    ov.push_back( Offset(0, 5) );
    ov.push_back( Offset(0) );
    ov.push_back( Offset(0, 3) );
    ov.push_back( Offset(0, 3, 6) );
    ov.push_back( Offset(1, 6, 3, 0) );
    ov.push_back( Offset(1, 6, 3) );
    ov.push_back( Offset(1, 6) );
    ov.push_back( Offset(1) );
    ov.push_back( Offset(1, 4) );
    ov.push_back( Offset(1, 4, 7) );
    ov.push_back( Offset(2, 7, 4, 1) );
    ov.push_back( Offset(2, 7, 4) );
    ov.push_back( Offset(2, 7) );
    ov.push_back( Offset(2) );
    ov.push_back( Offset(2, 5) );
    ov.push_back( Offset(2, 5, 0) );
    ov.push_back( Offset(3, 0, 5, 2) );
    ov.push_back( Offset(3, 0, 5) );
    ov.push_back( Offset(3, 0) );
    ov.push_back( Offset(3) );
    ov.push_back( Offset(3, 6) );
    ov.push_back( Offset(3, 6, 1) );
    ov.push_back( Offset(4, 1, 6, 3) );
    
    Offset r = ov[0];
    for (int i = 1; i < ov.size(); i++) {
        Offset oo = ov[i];
        int d = offset_distance(r, oo);
        
        cout << "r = " << r << ",\t oo = " << oo << ",\t distance = " << d << endl;
    }


    Offset o;
    cout << "r = " << r << ",\t oo = " << o << ",\t distance = " << offset_distance(r, o) << endl;
    cout << "r = " << o << ",\t oo = " << r << ",\t distance = " << offset_distance(o, r) << endl;
    cout << "r = " << o << ",\t oo = " << o << ",\t distance = " << offset_distance(o, o) << endl;

    /*
    o = Offset(0, 5, 2, 7);
    v = o.get_vector(); 
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(7, 2, 5, 0);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;

    o = Offset(0, 3, 6, 1);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(1, 6, 3, 0);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;

    o = Offset(1, 4, 7, 2);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(2, 7, 4, 1);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;


    o = Offset(2, 5, 0, 3);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(3, 0, 5, 2);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;


    o = Offset(3, 6, 1, 4);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(4, 1, 6, 3);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;


    o = Offset(4, 7, 2, 5);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(5, 2, 7, 4);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;


    o = Offset(5, 0, 3, 6);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(6, 3, 0, 5);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;


    o = Offset(6, 1, 4, 7);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    
    o = Offset(7, 4, 1, 6);
    v = o.get_vector();
    r.clear();
    Dehn_reduce_word(r, v);
    s = Offset(r);
    cout << o << " = " << s << endl;
    cout << "--------------" << endl;
    */

    return 0;
}

