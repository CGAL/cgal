#include <CGAL/box_intersection_d.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cassert>

struct Box {
    typedef int            NT;
    typedef std::ptrdiff_t ID;
    int lo[2], hi[2];
    Box( int lo0, int lo1, int hi0, int hi1) { lo[0]=lo0; lo[1]=lo1; hi[0]=hi0; hi[1]=hi1;}
    static int dimension() { return 2; }
    int min_coord(int dim) const { return lo[dim]; }
    int max_coord(int dim) const { return hi[dim]; }
    // id-function using address of current box,
    // requires to work with pointers to boxes later
    std::ptrdiff_t id() const { return (std::ptrdiff_t)(this); }
};

// 9 boxes of a grid
Box boxes[9] = { Box( 0,0,1,1),  Box( 1,0,2,1),  Box( 2,0,3,1), // low
                 Box( 0,1,1,2),  Box( 1,1,2,2),  Box( 2,1,3,2), // middle
                 Box( 0,2,1,3),  Box( 1,2,2,3),  Box( 2,2,3,3)};// upper
// 2 selected boxes as query; center and upper right
Box query[2] = { Box( 1,1,2,2),  Box( 2,2,3,3)};

// With the special id-function we need to work on box pointers
Box* b_ptr[9] = { boxes,   boxes+1, boxes+2, boxes+3, boxes+4, boxes+5,
                  boxes+6, boxes+7, boxes+8};
Box* q_ptr[2] = { query,   query+1};

// callback function object writing results to an output iterator
template <class OutputIterator>
struct Report {
    OutputIterator it;
    Report( OutputIterator i) : it(i) {} // store iterator in object
    // We write the position with respect to 'boxes' to the output iterator
    // assuming that box b (the query box) is not interesting in the result.
    void operator()( const Box* a, const Box*) {
        *it++ = ( reinterpret_cast<Box*>(a->id()) - boxes);
    }
};
template <class Iter> // helper function to create the function object
Report<Iter> report( Iter it) { return Report<Iter>(it); }

int main() {
    // run the intersection algorithm and store results in a vector
    std::vector<std::size_t> result;
    CGAL::box_intersection_d( b_ptr, b_ptr+9, q_ptr, q_ptr+2,
                              report( std::back_inserter( result)),
                              std::ptrdiff_t(0));
    // sort and check result
    std::sort( result.begin(), result.end());
    std::size_t chk[13] = {0,1,2,3,4,4,5,5,6,7,7,8,8};
    assert( result.size()==13 && std::equal(chk,chk+13,result.begin()));
    return 0;
}
