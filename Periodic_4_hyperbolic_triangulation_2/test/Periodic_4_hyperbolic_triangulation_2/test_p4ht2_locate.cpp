

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

using namespace CGAL;

typedef CORE::Expr                                              					NT;					
typedef CGAL::Cartesian<NT>                                     					Kernel;
typedef Kernel::FT                                                                  FT;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> 		Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>				Triangulation;
typedef Triangulation::Face_handle                                                  Face_handle;
typedef Triangulation::Locate_type                                                  Locate_type;
typedef Triangulation::Offset                                                       Offset;
typedef Triangulation::Point                                                        Point;

int ccw(int i) {
    return (i+1)%3;
}

std::ostream& operator<<(std::ostream& s, const Locate_type& lt) {
    switch(lt) {
        case Triangulation::VERTEX: s << "VERTEX";  break;
        case Triangulation::FACE:   s << "FACE";    break;
        case Triangulation::EDGE:   s << "EDGE";    break;
        default:                    s << "ZABABA";  break;
    }
    return s;
}


int main(void) {

    Triangulation tr;    
    tr.insert_dummy_points();

    assert(tr.is_valid());

    Locate_type lt;
    int li;
    Point query(0,0);

    Face_handle fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    FT F0 = FT(0);
    FT F1 = FT(1);
    FT F2 = FT(2);

    query = Point(FT( CGAL::sqrt(CGAL::sqrt(F2) - F1)), F0); //  μ(s_0)
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    query = Point(FT( CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2) ), FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)) );      //  μ(s_1)
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    query = Point(F0, FT( CGAL::sqrt(CGAL::sqrt(F2) - F1))); //  μ(s_2)
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    query = Point(-FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)), FT(CGAL::sqrt( (CGAL::sqrt(F2) - F1) / F2)) );       //  μ(s_3)
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    query = Point(FT( CGAL::sqrt(CGAL::sqrt(F2) + F1)/ F2 ), - FT( CGAL::sqrt(CGAL::sqrt(F2) - F1)/ F2) );           //  v_0
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);

    typedef Aff_transformation_2<Traits> Transformation; 
    Transformation rotate(ROTATION, FT(1)/CGAL::sqrt(FT(2)), FT(1)/CGAL::sqrt(FT(2))); 
    Transformation rot8(ROTATION, -CGAL::sqrt(FT(2)-CGAL::sqrt(FT(2)))/FT(2), CGAL::sqrt(FT(2)+CGAL::sqrt(FT(2)))/FT(2));
    Transformation shrink(SCALING, FT(FT(1)/FT(2)));

    FT i1(CGAL::sqrt(FT(2) - sqrt(FT(2))));
    FT i2(FT(2)*i1);
    FT i3(FT(2)*CGAL::sqrt(FT(2)) - i2);
    FT i4(-1 + i3);
    query = Point(FT( CGAL::sqrt( i4 )), FT(0));          // internal point

    query = rot8(query);
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::VERTEX);
    for (int i = 1; i < 8; i++) {
        query = rotate(query);
        fh = tr.locate(query, lt, li);
        cout << "Located " << query << " as " << lt << endl;
        assert(lt == Triangulation::VERTEX);
    }

    query = shrink(query);
    fh = tr.locate(query, lt, li);
    cout << "Located " << query << " as " << lt << endl;
    assert(lt == Triangulation::EDGE);
    for (int i = 1; i < 8; i++) {
        query = rotate(query);
        fh = tr.locate(query, lt, li);
        cout << "Located " << query << " as " << lt << endl;
        assert(lt == Triangulation::EDGE);
    }

    return 0;
}