#include <fstream>
#include <list>
#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_cartesian.h>
#include <CEP/Visibility_complex/Visibility_complex_segment_traits.h>
#include <CEP/Visibility_complex/Visibility_complex_items.h>
#include <CEP/Visibility_complex/Visibility_complex_2.h>

typedef CGAL::Gmpz   FT;
typedef CGAL::Simple_cartesian<FT>                      Rep;

typedef CGAL::Visibility_complex_segment_traits<Rep>    Gt;
typedef Gt::Disk                                        Segment;

// ---------------------------------------------------------------------

template <class Vc>
struct My_face : public CGAL::Visibility_complex_face_base<Vc>
{
    int size;
    My_face() : size(0) { }
};

struct My_items : public CGAL::Visibility_complex_items
{
    template <class Vc>
    struct Face_wrapper {
	typedef My_face<Vc>   Face;
    };
};

// ---------------------------------------------------------------------

typedef CGAL::Visibility_complex_2<Gt,My_items>         Visibility_complex;
typedef Visibility_complex::Antichain                   Antichain;
typedef Visibility_complex::Edge_handle                 Edge_handle;
typedef Visibility_complex::Vertex                      Vertex;

// ---------------------------------------------------------------------
// For a positive edge the three adjacent faces of e are
// dl(e) , ur(e) , ul(e)
// Otherwise the three adjacent faces are
// dl(e) , dr(e) , ul(e)

void increment(Edge_handle e)
{
    if (e->sign()) {
	++e->dl()->size;
	++e->ul()->size;
	++e->ur()->size;
    }
    else {
	++e->dl()->size;
	++e->dr()->size;
	++e->ul()->size;
    }
}

// ---------------------------------------------------------------------

int main()
{
    std::list<Segment> D;

    // Reading segments from file
    std::ifstream ifs("input");
    std::istream_iterator<Segment> ifs_it(ifs),ifs_end;
    std::copy(ifs_it,ifs_end,back_inserter(D));

    // Computing the initia; antichain to sweep the complex
    std::list<Vertex> V; // empty constraint list
    Antichain A(D.begin(),D.end(),V.begin(),V.end());

    Antichain::Linear_sweep_iterator v = A.sweep_begin();
    for ( ; v != A.sweep_end() ; ++v) 
    {
	// The two edges cw_source_edge(v) and cw_target(v) are being 
	// swept. We increment their three adjacent faces.
	increment(v->cw_source_edge());
	increment(v->cw_target_edge());

	// The face inf(v) has been completely swept. We print its size
	std::cout << v->inf()->size << std::endl;
    }
    return 0;
}
