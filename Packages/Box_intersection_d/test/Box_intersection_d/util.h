#ifndef BOX_INTERSECTION_D_UTIL_H
#define BOX_INTERSECTION_D_UTIL_H

#define USING(T) typedef typename Defs::T T;
    USING(Traits)
    USING(Box)
    USING(Box_container)
    USING(Box_pair)
    USING(Result_container)
#undef USING

static void readBoxesFromFile( FILE *infile, Box_container& boxes ) {
  int numBoxes, numDim;
  int boxNum, dim;

  fscanf(infile, "%d %d\n", &numBoxes, &numDim);
  std::vector< int > min( numDim ), max( numDim );
  /* Read boxes */
  for(boxNum = 0; boxNum < numBoxes; boxNum++) {
      for(dim = 0; dim < numDim; dim++)
          fscanf( infile, "[%d, %d) ", &min[dim], &max[dim] );
      boxes.push_back( Box( &min[0], &max[0] ) );
      fscanf(infile, "\n");
  }
}

static void fill_boxes( unsigned int n, Box_container& boxes ) {
    NT maxEdgeLength = (NT) pow(n, (DIM-1.0)/DIM);

    for( unsigned int i = 0; i < n; ++i ) {
        NT lo[DIM], max[DIM];
        for( unsigned int d = 0; d < DIM; ++d ) {
            lo[d] =
                (NT)(CGAL::default_random.get_double() * (n - maxEdgeLength));
            max[d] =
        (NT)(lo[d] + 1 + (CGAL::default_random.get_double() * maxEdgeLength));
        }
        boxes.push_back( Box( &lo[0], &max[0]) );
    }
}

static void assert_intersection( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < DIM; ++dim ) {
        if( Traits::does_intersect( a, b, dim ) == false ) {
            std::cout << "does not intersect!" << std::endl;
            //cout << a << endl << b << endl;
            exit(-1);
        }
    }
}

struct Counter_callback {
    unsigned int counter;
    Counter_callback() : counter( 0 ) {}
    void operator()( const Box& a, const Box& b ) {
        assert_intersection( a, b );
        ++counter;
    }
};

template< class Storage >
struct Storage_callback : public Counter_callback {
    Storage& storage;
    Storage_callback( Storage& storage ) : storage( storage ) {}
    void operator()( const Box& a, const Box& b ) {
        Counter_callback::operator()(a,b);
        storage.push_back( std::make_pair( a, b ) );
    }
};

static bool
areEqual( const Box& a, const Box& b ) {
    for( unsigned int dim = 0; dim < DIM; ++dim )
        if( Traits::min( a, dim ) != Traits::min( b, dim ) ||
            Traits::max( a, dim ) != Traits::max( b, dim )   )
            return false;
    return true;
}

static bool
areEqual( const Box_pair& a, const Box_pair& b ) {
    return( areEqual(a.first,b.first)  && areEqual(a.second,b.second) ||
            areEqual(a.first,b.second) && areEqual(a.second,b.first ));
}

template< class Storage >
static unsigned int countMissingItems( Storage& a, Storage& b ) {
    unsigned int missing = 0;
    for( typename Storage::iterator it = a.begin(); it != a.end(); ++it )
        for( typename Storage::iterator it2 = b.begin(); it2 != b.end(); ++it2 )
            if( areEqual(*it,*it2) )
                ++missing;
    return missing;
}

template< class Storage >
static unsigned int countDuplicates( Storage& storage ) {
    unsigned int counter = 0;
    typedef typename Storage::iterator IT;
    for( IT it = storage.begin(); it != storage.end(); ++it )
        for( IT it2 = it; it2 != storage.end(); ++it2 )
            if( it != it2 &&  areEqual(*it,*it2) ) {
                //cout << it->first.id() << " <-> "
                //     << it->second.id() << endl;
                ++counter;
            }
    return counter;
}

#endif
