#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_2.h>
#include <CGAL/Visibility_complex_segment_traits.h>

typedef CGAL::Simple_cartesian<int>                  Rep;
typedef CGAL::Visibility_complex_segment_traits<Rep> Gt;
typedef Gt::Disk                                     Segment;
typedef CGAL::Visibility_complex_2<Gt>               Visibility_complex;
typedef Visibility_complex::Antichain                Antichain;
typedef Visibility_complex::Linear_sweep_iterator    Linear_sweep_iterator;

int main()
{
    // Reading circles from file
    std::list<Segment> D;
    std::ifstream ifs("input");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,back_inserter(D));

    // Computing the antichain
    Antichain A(D.begin(),D.end());

    // Sweep of the visibility complex using an iterator
    Linear_sweep_iterator v(&A), vend(&A,0);

    for ( ; v != vend ; ++v) 
	std::cout << *v << std::endl;

    return 0;
}
