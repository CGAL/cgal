// enable invariant checking
#define SEGMENT_TREE_CHECK_INVARIANTS 1
#include <CGAL/Box_intersection_d.h>

#include <CGAL/Box_intersection_d/all_pairs.h>
#include <CGAL/Box_intersection_d/one_way_scan.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <cstdio>
#include <cmath>

bool test_failed = false;


typedef CGAL::Box_intersection_d::Box_d< int, 3 >     Box;
typedef CGAL::Box_intersection_d::Box_traits_d< Box > Box_adapter;
typedef CGAL::Box_intersection_d::Box_predicate_traits_d<
                                                  Box_adapter, true > Traits;
typedef std::vector< Box >      Box_container;
typedef std::pair< Box, Box >   Box_pair;
typedef std::vector< Box_pair > Result_container;


static void readBoxesFromFile( FILE *infile, Box_container& boxes )
{
  int numBoxes, numDim;
  int boxNum, dim;

  fscanf(infile, "%d %d\n", &numBoxes, &numDim);
  std::vector< int > lo( numDim ), hi( numDim );
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
            std::cout << "does not intersect!" << std::endl;
            //cout << a << endl << b << endl;
            exit(-1);
        }
    }
}

template< class Storage >
struct Storage_callback {
    unsigned int counter;
    Storage& storage;
    Storage_callback( Storage& storage ) : counter( 0 ), storage( storage ) {}
    void operator()( const Box& a, const Box& b ) {
        assertIntersection( a, b );
        ++counter;
        storage.push_back( std::make_pair( a, b ) );
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
operator==( const Box_pair& a, const Box_pair& b ) {
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
                //cout << it->first.id() << " <-> "
                //     << it->second.id() << endl;
                ++counter;
            }
    return counter;
}

static void
test( const char* filename1, const char* filename2 )
{
    Box_container boxes1, boxes2;
    Result_container result_all_pairs, result_tree;
    FILE *infile1, *infile2;
    infile1 = fopen( filename1, "r");
    infile2 = fopen( filename2, "r");

    readBoxesFromFile( infile1, boxes1 );
    readBoxesFromFile( infile2, boxes2 );

    std::cout << std::endl;
    Storage_callback< Result_container >
        callback1( result_all_pairs ),
        callback2( result_tree );

    std::cout << "all pairs ...... " << std::flush;
    CGAL::Timer timer;
    timer.start();
    CGAL::Box_intersection_d::all_pairs( boxes1.begin(), boxes1.end(),
                                         boxes2.begin(), boxes2.end(),
                                         callback1, Traits(), 2 );
    timer.stop();
    std::cout << "got " << callback1.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;


    std::cout << "one way scan ...... " << std::flush;
    timer.reset();
    timer.start();
    CGAL::Box_intersection_d::one_way_scan( boxes1.begin(), boxes1.end(),
                                            boxes2.begin(), boxes2.end(),
                                            callback2, Traits(), 2 );
    CGAL::Box_intersection_d::one_way_scan( boxes2.begin(), boxes2.end(),
                                            boxes1.begin(), boxes1.end(),
                                            callback2, Traits(), 2 );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;
    callback2.counter = 0;
    result_tree.clear();

    std::cout << "segment tree ... " << std::flush;
    timer.reset();
    timer.start();
    const unsigned int n = boxes1.size();
    const unsigned int cutoff = n < 2000 ? 6 : n / 100;
    CGAL::box_intersection_d_custom( boxes1.begin(), boxes1.end(),
                                     boxes2.begin(), boxes2.end(),
                                     callback2, Traits(), cutoff );
    timer.stop();
    std::cout << "got " << callback2.counter << " intersections in "
              << timer.time() << " seconds." << std::endl;

    if( callback1.counter != callback2.counter ) {
        unsigned int missing    = countMissingItems( result_all_pairs,
                                                     result_tree );
        unsigned int duplicates = countDuplicates( result_tree );
        std::cout << "!! failed !! " << missing  << " missing and "
             << duplicates << " duplicate intersections in tree result."
             << std::endl;
        test_failed = true;
    }
    else
        std::cout << "--- passed --- " << std::endl;
    fclose( infile1 );
    fclose( infile2 );
}


int main( int argc, char ** argv ) {
    for( unsigned int n = 1; n <= 6; ++n ) {
        std::stringstream file1, file2;
        file1 << "data/test" << n << "_set1.box" << std::ends;
        file2 << "data/test" << n << "_set2.box" << std::ends;
        test( file1.str().c_str(), file2.str().c_str() );
    }
    if ( test_failed)
        return 1;
    return 0;
}

