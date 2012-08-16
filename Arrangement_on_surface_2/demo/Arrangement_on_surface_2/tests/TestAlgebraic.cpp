#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

typedef CGAL::CORE_algebraic_number_traits Nt_traits;
typedef Nt_traits::Integer Coefficient;
typedef CGAL::Arr_algebraic_segment_traits_2< Coefficient > Traits;
typedef Traits::Algebraic_kernel_d_1 Algebraic_kernel_d_1;
typedef Traits::Algebraic_kernel_d_2 Algebraic_kernel_d_2;
typedef Traits::Algebraic_real_1 Algebraic_real_1;
typedef Traits::Construct_curve_2 Construct_curve_2;
typedef Traits::Curve_2 Curve_2;
typedef Traits::X_monotone_curve_2 X_monotone_curve_2;
typedef Traits::Polynomial_2 Polynomial_2;
typedef Traits::Bound Bound;
typedef Traits::Point_2 Point_2;
typedef Traits::Multiplicity Multiplicity;
typedef Traits::Make_x_monotone_2 Make_x_monotone_2;
typedef Traits::Intersect_2 Intersect_2;
//typedef Algebraic_kernel_d_2::Make

typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_kernel_d_2 > CKvA_2;

Polynomial_2 makeParabola( )
{
    Polynomial_2 x = CGAL::shift( Polynomial_2( 1 ), 1, 0 );
    Polynomial_2 y = CGAL::shift( Polynomial_2( 1 ), 1, 1 );
    Polynomial_2 parabola = y - x*x;

    return parabola;
}

X_monotone_curve_2 makeVerticalLine( Bound x )
{
    Traits traits;
    Traits::Construct_point_2 constructPoint =
        traits.construct_point_2_object( );
    Traits::Construct_x_monotone_segment_2 constructSegment = 
        traits.construct_x_monotone_segment_2_object( );

    std::vector< X_monotone_curve_2 > curves;
    Point_2 p1 = constructPoint( Algebraic_real_1(x), Algebraic_real_1(Bound( -10000 )) );
    Point_2 p2 = constructPoint( x, Bound( +10000 ) );
    constructSegment( p1, p2, std::back_inserter( curves ) );
    return curves[ 0 ];
}

typedef CGAL::Cartesian< Coefficient > Kernel;
typedef Kernel::Point_2 Kernel_point_2;

int main( )
{
    Algebraic_real_1 real( 1 );
    //CGAL::Qt::Converter< Algebraic_kernel_d_2 > testConverter;
    //CGAL::Qt::Converter< CKvA_2 > testConverter;
    
    //CGAL::Qt::Converter< Cartesian > testConverter;
    Kernel_point_2 testPt( 1, 2 );
    Point_2 testPt2( testPt.x( ), testPt.y( ) );
    Traits traits;
    Construct_curve_2 constructCurve = traits.construct_curve_2_object( );
    Curve_2 curve = constructCurve( makeParabola( ) );
    Make_x_monotone_2 mm = traits.make_x_monotone_2_object( );
    std::vector< CGAL::Object > curves;
    mm( curve, std::back_inserter( curves ) );
    std::cout << curves.size( ) << std::endl;
    X_monotone_curve_2 c1;
    CGAL::assign( c1, curves[ 0 ] );
    double lb = -3;
    double ub = 3;
    double step = 6.0 / 1000;

    for ( int i = 0; i < 1000; ++i )
    {
        X_monotone_curve_2 c2 = makeVerticalLine( lb + step * i );

        CGAL::Object o;
        CGAL::Oneset_iterator< CGAL::Object > oi( o );
        Intersect_2 intersect = traits.intersect_2_object( );
        intersect( c1, c2, oi );
        std::pair< Point_2, Multiplicity > res;
        CGAL::assign( res, o );
        std::pair< double, double > approx = res.first.to_double( );
        std::cout << approx.first << " " << approx.second << std::endl;
    }


    return 0;
}
