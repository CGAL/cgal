#ifndef CGAL_BOX_INTERSECTION_D_UTIL_H
#define CGAL_BOX_INTERSECTION_D_UTIL_H

#include <vector>
#include <algorithm> // for pair
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <CGAL/Random.h>

template< class NT, int DIM, bool CLOSED = true >
struct Util {
    typedef NT Number_type;
    typedef CGAL::Box_intersection_d::Box_d< Number_type, DIM >  Box;
    typedef CGAL::Box_intersection_d::Box_traits_d< Box >        B_traits;
    typedef CGAL::Box_intersection_d::Predicate_traits_d< B_traits, CLOSED > 
                                                                 Traits;
    typedef std::vector< Box >      Box_container;
    typedef std::pair< Box, Box >   Box_pair;
    typedef std::vector< Box_pair > Result_container;

    static void readBoxesFromFile( std::FILE *infile, Box_container& boxes ) {
      int numBoxes, numDim;
      int boxNum, dim;
    
      std::fscanf(infile, "%d %d\n", &numBoxes, &numDim);
      std::vector< int > minc( numDim ), maxc( numDim );
      /* Read boxes */
      for(boxNum = 0; boxNum < numBoxes; boxNum++) {
          for(dim = 0; dim < numDim; dim++)
              std::fscanf( infile, "[%d, %d) ", &minc[dim], &maxc[dim] );
          boxes.push_back( Box( &minc[0], &maxc[0] ) );
          std::fscanf(infile, "\n");
      }
    }
    
    static void fill_boxes( unsigned int n, Box_container& boxes ) {
        NT maxEdgeLength = (NT) std::pow(n, (DIM-1.0)/DIM);
    
        for( unsigned int i = 0; i < n; ++i ) {
            NT lo[DIM], max[DIM];
            for( int d = 0; d < DIM; ++d ) {
                lo[d] =
                    (NT)(CGAL::get_default_random().get_double() * (n - maxEdgeLength));
                max[d] =
            (NT)(lo[d] + 1 + (CGAL::get_default_random().get_double() * maxEdgeLength));
            }
            boxes.push_back( Box( &lo[0], &max[0]) );
        }
    }
    
    static void assert_intersection( const Box& a, const Box& b ) {
        for( int dim = 0; dim < DIM; ++dim ) {
            if( Traits::does_intersect( a, b, dim ) == false ) {
                std::cout << "does not intersect!" << std::endl;
                //cout << a << endl << b << endl;
                std::exit(-1);
            }
        }
    }
    
    struct Counter_callback {
        unsigned int& counter;
        Counter_callback(unsigned int& i)
          : counter(i)
      {}
        void operator()( const Box& a, const Box& b ) {
            assert_intersection( a, b );
            ++counter;
        }
        unsigned int get_counter() { return counter; }
        void reset_counter() { counter = 0; }
    };
    
    template< class Storage = Result_container >
    struct Storage_callback : public Counter_callback {
        Storage& storage;
      Storage_callback( Storage& storage, unsigned int& i ) : Counter_callback(i), storage( storage ) {}
        void operator()( const Box& a, const Box& b ) {
            Counter_callback::operator()(a,b);
            storage.push_back( std::make_pair( a, b ) );
        }
    };
    
    static bool
    areEqual( const Box& a, const Box& b ) {
        for( int dim = 0; dim < DIM; ++dim )
            if( Traits::min_coord( a, dim ) != Traits::min_coord( b, dim ) ||
                Traits::max_coord( a, dim ) != Traits::max_coord( b, dim )   )
                return false;
        return true;
    }
    
    static bool
    areEqual( const Box_pair& a, const Box_pair& b ) {
        return ( areEqual(a.first,b.first)  && areEqual(a.second,b.second) ) ||
               ( areEqual(a.first,b.second) && areEqual(a.second,b.first ) );
    }

    template< class Storage >
    static unsigned int
    countMissingItems( Storage& a, Storage& b ) {
        unsigned int missing = 0;
        typedef typename Storage::iterator iterator;
        for( iterator it = a.begin(); it != a.end(); ++it )
            for( iterator it2 = b.begin(); it2 != b.end(); ++it2 )
                if( areEqual(*it,*it2) )
                    ++missing;
        return missing;
    }
    
    template< class Storage >
    static unsigned int
    countDuplicates( Storage& storage ) {
        unsigned int counter = 0;
        typedef typename Storage::iterator iterator;
        for( iterator it = storage.begin(); it != storage.end(); ++it )
            for( iterator it2 = it; it2 != storage.end(); ++it2 )
                if( it != it2 &&  areEqual(*it,*it2) ) {
                    //cout << it->first.id() << " <-> "
                    //     << it->second.id() << endl;
                    ++counter;
                }
        return counter;
    }


};


#endif
