#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_segment_traits.h>
#include <CGAL/Visibility_complex_items.h>
#include <CGAL/Visibility_complex_2.h>

typedef CGAL::Gmpz   FT;
typedef CGAL::Simple_cartesian<FT>                      Rep;

typedef CGAL::Visibility_complex_segment_traits<Rep>    Gt;
typedef Gt::Disk                                        Segment;

// ---------------------------------------------------------------------

template <class V>
struct My_face : public CGAL::Visibility_complex_face_base<V>
{
    int size;
    My_face() : size(0) { }
};

struct My_items : public CGAL::Visibility_complex_items
{
    template <class V>
    struct Face_wrapper {
	typedef My_face<V>   Face;
    };
};

// ---------------------------------------------------------------------

typedef CGAL::Visibility_complex_2<Gt,My_Items>         Visibility_complex;
typedef Visibility_complex::Antichain                   Antichain;
typedef Visibility_complex::Linear_sweep_iterator       Linear_sweep_iterator;
typedef Visibility_complex::Edge_handle                 Edge_handle;

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
    Antichain A(D.begin(),D.end());
    Linear_sweep_iterator v(&A) , vend(&A,0);

    for ( ; v != vend ; ++v) 
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
