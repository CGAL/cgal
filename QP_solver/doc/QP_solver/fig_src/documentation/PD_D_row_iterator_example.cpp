template < class Point >
class PD_D_row_iterator {

  // ...

  // data members
  Point_it      points;       // iterator referring to C matrix
  int           i;            // row index
  int           j;            // column index

public:
  PD_D_row_iterator( Point_it it, int row) : points(it), i(row), j(0) { }

  CT  operator * ( ) {
    return std::inner_product( points[ i].begin( ), points[ i].end( ),
      points[ j].begin( ), CT( ));
  }

  // other operators required for random access iterators
  // ...

};

