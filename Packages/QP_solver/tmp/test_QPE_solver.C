// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : test/_QP_solver/test_QP_solver.C
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/21
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for the QP solver
// ============================================================================

#include <CGAL/QPE_solver.h>
#include <CGAL/QPE_full_exact_pricing.h>
#include <CGAL/QPE_partial_exact_pricing.h>
#include <CGAL/QPE_full_filtered_pricing.h>
#include <CGAL/QPE_partial_filtered_pricing.h>
#include <CGAL/_QP_solver/Double.h>

#include <functional>
#include <iostream>

typedef  std::vector<double>  Vector;
typedef  std::vector<Vector>  Matrix;

struct Begin
    : public CGAL_STD::unary_function< Vector, Vector::const_iterator > {
    result_type  operator () ( const Vector& v) const { return v.begin(); }
};

typedef  CGAL::QPE_transform_iterator_1< Matrix::const_iterator, Begin>
                                        Vector_iterator;
typedef  Vector::const_iterator         Entry_iterator;

struct Rep {
    typedef  GMP::Double  ET;

    typedef  Vector_iterator  A_iterator;
    typedef   Entry_iterator  B_iterator;
    typedef   Entry_iterator  C_iterator;
    typedef  Vector_iterator  D_iterator;

    enum Row_type { LESS_EQUAL = -1, EQUAL, GREATER_EQUAL};
    typedef  Row_type*  Row_type_iterator;
    //    typedef  CGAL::QPE_const_value_iterator<Row_type>  Row_type_iterator;
    

    typedef  CGAL::Tag_false  Is_linear;
    typedef  CGAL::Tag_true  Is_symmetric;
    typedef  CGAL::Tag_true  Has_no_inequalities;
};

int
main( int argc, char** argv)
{
    int verbose = 0;
    if ( argc > 1) verbose = atoi( argv[ 1]);

//     {
//       Matrix  A( 3);
//       Vector  b;
//       Vector  c;
//       Matrix  D( 3);

//       A[ 0].push_back( 3); A[ 1].push_back( 2); A[ 2].push_back( 0);
//       A[ 0].push_back( 0); A[ 1].push_back( 2); A[ 2].push_back( 1);

//       b.push_back( 1);
//       b.push_back( 2);

//       c.push_back( 1); c.push_back( 2); c.push_back( 3);
//       /*
// 	D[ 0].push_back( 3); D[ 0].push_back( 1); D[ 0].push_back( -9);
// 	D[ 1].push_back( 1); D[ 1].push_back( -2); D[ 1].push_back( 2);
// 	D[ 2].push_back( -9); D[ 2].push_back( 2); D[ 2].push_back( 1);
//       */
//       //*
//       D[ 0].push_back(  2); D[ 0].push_back( -5); D[ 0].push_back( -5);
//       D[ 1].push_back( -5); D[ 1].push_back( 11); D[ 1].push_back( -12);
//       D[ 2].push_back( -5); D[ 2].push_back( -12); D[ 2].push_back( -11);
//       //*/

//       Rep::Row_type*  row_types = new Rep::Row_type[ 2];
//       row_types[ 0] = Rep::EQUAL;
//       row_types[ 1] = Rep::EQUAL;

//       CGAL::QPE_solver<Rep>              solver;
//       solver.set_verbosity( verbose);
//       solver.set( 3, 2,
// 		  Vector_iterator( A.begin(), Begin()), b.begin(),
// 		  c.begin(), Vector_iterator( D.begin(), Begin()),
// 		  row_types);
//       {
// 	CGAL::QPE_full_exact_pricing<Rep>  strategy;
// 	solver.set_pricing_strategy( strategy);
// 	solver.init();
// 	solver.solve();
// 	std::cerr << "-----------------------------------------------------\n";
//       }
//       {
// 	CGAL::QPE_partial_exact_pricing<Rep>  strategy;
// 	solver.set_pricing_strategy( strategy);
// 	solver.init();
// 	solver.solve();
// 	std::cerr << "-----------------------------------------------------\n";
//       }
//       {
// 	CGAL::QPE_full_filtered_pricing<Rep>  strategy;
// 	solver.set_pricing_strategy( strategy);
// 	solver.init();
// 	solver.solve();
// 	std::cerr << "-----------------------------------------------------\n";
//       }
//       {
// 	CGAL::QPE_partial_filtered_pricing<Rep>  strategy;
// 	solver.set_pricing_strategy( strategy);
// 	solver.init();
// 	solver.solve();
// 	std::cerr << "-----------------------------------------------------\n";
//       }
//   }

    {
      Matrix  A( 2);
      Vector  b;
      Vector  c;
      Matrix  D( 2);

      A[ 0].push_back( 1); A[ 1].push_back( -1); 

      b.push_back( 0);
 
      c.push_back( 0); c.push_back( 0);
      /*
	D[ 0].push_back( 3); D[ 0].push_back( 1); D[ 0].push_back( -9);
	D[ 1].push_back( 1); D[ 1].push_back( -2); D[ 1].push_back( 2);
	D[ 2].push_back( -9); D[ 2].push_back( 2); D[ 2].push_back( 1);
      */
      //*
      D[ 0].push_back(  1); D[ 0].push_back( -1); 
      D[ 1].push_back( -1); D[ 1].push_back(  1); 
      //*/

      Rep::Row_type*  row_types = new Rep::Row_type[ 1];
      row_types[ 0] = Rep::EQUAL;

      CGAL::QPE_solver<Rep>              solver;
      solver.set_verbosity( verbose);
      solver.set( 2, 1,
		  Vector_iterator( A.begin(), Begin()), b.begin(),
		  c.begin(), Vector_iterator( D.begin(), Begin()),
		  row_types);
      {
	CGAL::QPE_full_exact_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_partial_exact_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_full_filtered_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
      {
	CGAL::QPE_partial_filtered_pricing<Rep>  strategy;
	solver.set_pricing_strategy( strategy);
	solver.init();
	solver.solve();
	std::cerr << "-----------------------------------------------------\n";
      }
    }

    return 0;
}

// ===== EOF ==================================================================
