#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>
#include <CEP/Visibility_complex/Visibility_complex_segment_traits.h>

typedef CGAL::Simple_cartesian<int>                  Rep;
typedef CGAL::Visibility_complex_segment_traits<Rep> Gt;
typedef Gt::Disk                                     Segment;
typedef CGAL::Visibility_complex_2<Gt>               Visibility_complex;
typedef Visibility_complex::Antichain                Antichain;
typedef Antichain::Linear_sweep_iterator             Linear_sweep_iterator;
typedef Antichain::Vertex Vertex;

int main()
{
    // Reading circles from file
    std::list<Segment> D;
    std::ifstream ifs("input");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,back_inserter(D));

    // Computing the antichain
    std::list<Vertex> V; // dummy empty constraint list
    Antichain A(D.begin(),D.end(),V.begin(),V.end());

    // Sweep of the visibility complex using an iterator
    Linear_sweep_iterator v = A.sweep_begin();

    for ( ; v != A.sweep_end() ; ++v) 
	std::cout << *v << std::endl;

    return 0;
}
