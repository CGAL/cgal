//#include "Box_intersection_d.h"
#include <CGAL/Box_intersection_d/box_traits.h>

#include <CGAL/Box_intersection_d/all_pairs.h>
#include <CGAL/Box_intersection_d/one_way_scan.h>

// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d/segment_tree.h>

#include "Timer.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>


using namespace std;

typedef CGAL::Default_Bbox_d< int, 3 >         Box;
typedef CGAL::Default_Bbox_d_Adapter< Box >    BoxAdapter;
typedef CGAL::Default_Box_Traits< BoxAdapter > Traits;
typedef vector< Box >     BoxContainer;
typedef pair< Box, Box >  BoxPair;
typedef vector< BoxPair > ResultContainer;


static void readBoxesFromFile( FILE *infile, BoxContainer& boxes )
{
  int numBoxes, numDim;
  int boxNum, dim;

  fscanf(infile, "%d %d\n", &numBoxes, &numDim);
  vector< int > lo( numDim ), hi( numDim );
  /* Read boxes */
  for(boxNum = 0; boxNum < numBoxes; boxNum++) {
      for(dim = 0; dim < numDim; dim++)
          fscanf( infile, "[%d, %d) ", &lo[dim], &hi[dim] );
      boxes.push_back( Box( &lo[0], &hi[0] ) );
      fscanf(infile, "\n");
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
    for( typename Storage::iterator it = a.begin(); it != a.end(); ++it )
        for( typename Storage::iterator it2 = b.begin(); it2 != b.end(); ++it2 )
            if( *it == *it2 )
                ++missing;
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
test( const char* filename1, const char* filename2 )
{
    BoxContainer boxes1, boxes2;
    ResultContainer result_all_pairs, result_tree;
    FILE *infile1, *infile2;
    infile1 = fopen( filename1, "r");
    infile2 = fopen( filename2, "r");

    readBoxesFromFile( infile1, boxes1 );
    readBoxesFromFile( infile2, boxes2 );

    cout << endl;
    StorageCallback< ResultContainer >
        callback1( result_all_pairs ),
        callback2( result_tree );

    cout << "all pairs ...... " << flush;
    Timer timer;
    timer.start();
    CGAL::all_pairs( boxes1.begin(), boxes1.end(),
                     boxes2.begin(), boxes2.end(), callback1, Traits(), 2 );
    timer.stop();
    cout << "got " << callback1.counter << " intersections in "
         << timer.t << " seconds." << endl;

    cout << "segment tree ... " << flush;
    timer.reset();
    timer.start();
    unsigned int n = boxes1.size();
    Traits::cutoff = n < 2000 ? 6 : n / 100;
    CGAL::segment_tree( boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(), callback2, Traits() );
    timer.stop();
    cout << "got " << callback2.counter << " intersections in "
         << timer.t << " seconds." <<endl;

    if( callback1.counter != callback2.counter ) {
        unsigned int missing    = countMissingItems( result_all_pairs,
                                                     result_tree );
        unsigned int duplicates = countDuplicates( result_tree );
        cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result."
             << endl;
    }
    else
        cout << "--- passed --- " << endl;
    fclose( infile1 );
    fclose( infile2 );
}


int main( int argc, char ** argv ) {
    for( unsigned int n = 1; n <= 6; ++n ) {
        stringstream file1, file2;
        file1 << "data/test" << n << "_set1.box" << ends;
        file2 << "data/test" << n << "_set2.box" << ends;
        test( file1.str().c_str(), file2.str().c_str() );
    }
}

