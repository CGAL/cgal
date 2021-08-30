template < class Point >
class PD_D_iterator {

  // ...

  // data members
  Point_it      points;       // iterator referring to C matrix
  int           i;            // row index

public:
  PD_D_iterator( Point_it it) : points( it), i( 0) { }

  Row_it  operator * ( ) { return Row_it( points, i);}

  // other operators required for random access iterators
  // ...
};

typedef  CGAL::QP_transform_iterator_1< Matrix::const_iterator>
                                                 Vector_iterator;

struct QPSolverTraits {
    typedef  GMP::Double                         ET;
    typedef  Vector_iterator                     A_iterator;
    typedef  CGAL::QP_const_value_iterator<IT>  B_iterator;
    typedef  CGAL::QP_const_value_iterator<IT>  C_iterator;
    typedef  PD_D_iterator<Point>                D_iterator;

    enum Row_type { LESS_EQUAL, EQUAL, GREATER_EQUAL};
    typedef  CGAL::QP_const_value_iterator<Row_type>  Row_type_iterator;

    typedef  CGAL::Tag_false  Is_linear;
    typedef  CGAL::Tag_true   Is_symmetric;
    typedef  CGAL::Tag_true   Has_full_row_rank;
    typedef  CGAL::Tag_false  Use_perturbation;
};
