#include <CGAL/Cartesian.h>
#include <CGAL/Visibility_complex_2.h>
#include <CGAL/Visibility_complex_circle_traits.h>

typedef CGAL::Cartesian<double>                     Rep;
typedef CGAL::Visibility_complex_circle_traits<Rep> Traits;
typedef Traits::Disk                                Circle;
typedef CGAL::Visibility_complex<Traits>            VComplex;

int main()
{
    std::list<Circle> D;

    // Reading circles from file
    std::ifstream ifs("input");
    std::istream_iterator<Circle> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,back_inserter(D));

    // Computing Visibibility Complex
    VComplex V(D.begin(),D.end());

    // Sice of V
    int size = 0;
    distance(V.vertices_begin(),V.vertices_end(),size);
    std::cout << size << std::endl;

    return 0;
}
