#include "bbox.h"
#include <bbox/box_traits.h>

#include <bbox/one_way_scan.h>

// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <bbox/segment_tree.h>

#include "Timer.h"

#include <iostream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>


using namespace std;

typedef double NumberType;
typedef Bbox_3< NumberType > Box;
typedef Bbox_3_Adapter< Box > BoxAdapter;
typedef Default_Box_Traits< BoxAdapter > Traits;
typedef vector< Box > BoxContainer;
typedef pair< Box, Box > BoxPair;
typedef vector< BoxPair > ResultContainer;


static void fill_boxes( unsigned int n, BoxContainer& boxes ) {
    unsigned int maxEdgeLength = (int) pow(n, (double)2/3);

    for( unsigned int i = 0; i < n; ++i ) {
        NumberType lo[3], hi[3];
        for( unsigned int d = 0; d < 3; ++d ) {
            lo[d] = (NumberType)(drand48() * (n - maxEdgeLength));
            hi[d] = lo[d] + 1 + (NumberType)(drand48() * maxEdgeLength);
        }
        boxes.push_back( Box( lo[0], lo[1], lo[2], hi[0], hi[1], hi[2] ) );
    }
}

static void assertIntersection( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < 3; ++dim ) {
        if( Traits::does_intersect( a, b, dim ) == false ) {
            cout << "does not intersect!" << endl;
            //cout << a << endl << b << endl;
            exit(-1);
        }
    }
}


template< class Storage >
struct StorageCallback {
    unsigned int counter;
    Storage& storage;
    StorageCallback( Storage& storage ) : counter( 0 ), storage( storage ) {}
    void operator()( const Box& a, const Box& b ) {
        assertIntersection( a, b );
        ++counter;
        storage.push_back( make_pair( a, b ) );
    }
};




bool
operator==( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < 3; ++dim )
        if( Traits::get_lo( a, dim ) != Traits::get_lo( b, dim ) ||
            Traits::get_hi( a, dim ) != Traits::get_hi( b, dim )   )
            return false;
    return true;
}

bool
operator==( const BoxPair& a, const BoxPair& b ) {
    return( a.first == b.first && a.second == b.second ||
            a.first == b.second && a.second == b.first );
}

template< class Storage >
unsigned int countMissingItems( Storage& a, Storage& b ) {
    unsigned int missing = 0;
    for( typename Storage::iterator it = a.begin(); it != a.end(); ++it ) {
        if( find( b.begin(), b.end(), *it ) == b.end() ) {
            ++missing;
            //cout << it->first << it->second << endl;
        }
    }
    return missing;
}

template< class Storage >
unsigned int countDuplicates( Storage& storage ) {
    unsigned int counter = 0;
    typedef typename Storage::iterator IT;
    for( IT it = storage.begin(); it != storage.end(); ++it )
        for( IT it2 = it; it2 != storage.end(); ++it2 )
            if( it != it2 &&  *it == *it2 ) {
                //cout << it->first.num() << " <-> "
                //     << it->second.num() << endl;
                ++counter;
            }
    return counter;
}

static void
test_n( unsigned int n )
{
    BoxContainer boxes1, boxes2;
    ResultContainer result_scanner, result_tree;
    cout << "generating random box sets with size " << n << " ... " << flush;
    fill_boxes( n, boxes1 );
    fill_boxes( n, boxes2 );
    cout << endl;
    StorageCallback< ResultContainer >
        callback1( result_scanner ),
        callback2( result_tree );

    cout << "one way scan ... " << flush;
    Timer timer;
    timer.start();
    one_way_scan( boxes1.begin(), boxes1.end(),
                  boxes2.begin(), boxes2.end(), callback1, Traits(), 2 );
    one_way_scan( boxes2.begin(), boxes2.end(),
                  boxes1.begin(), boxes1.end(), callback1, Traits(), 2 );
    timer.stop();
    cout << "got " << callback1.counter << " intersections in "
         << timer.t << " seconds."
         << endl;

    cout << "segment tree ... " << flush;
    timer.reset();
    timer.start();
    Traits::cutoff = n < 2000 ? 6 : n / 100;
    //Traits::cutoff = 5;
    segment_tree( boxes1.begin(), boxes1.end(),
                  boxes2.begin(), boxes2.end(), callback2, Traits(), 2 );
    timer.stop();
    cout << "got " << callback2.counter << " intersections in "
         << timer.t << " seconds." <<endl;

    if( callback1.counter != callback2.counter ) {
        unsigned int missing    = countMissingItems( result_scanner,
                                                     result_tree );
        unsigned int duplicates = countDuplicates( result_tree );
        cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result."
             << endl;
    }
    else
        cout << "--- passed --- " << endl;
}


int main( int argc, char ** argv ) {
    for( unsigned int n = 8; n < 500000; n = (int)(n * 1.3))
        test_n( n );
}

