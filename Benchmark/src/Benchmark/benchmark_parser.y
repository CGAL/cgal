/**************************************************************************
// Copyright (c) 2004  Max-Planck-Institut Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s) : Lutz Kettner
//             Franziska Ebert <febert@mpi-sb.mpg.de>
**************************************************************************/

%{
/* C/C++ declaration section */
/* ========================= */
#include <cassert>
#include <cstdlib>  /* for atoi */
#include <string>    /* for std::string */
#include <fstream>   /* for std::ifstream */
#include <sstream>   // for better error messages
#include <iostream>
#include <stack>     // for type checking

#include <CGAL/Benchmark/config.hpp>
#include <CGAL/Benchmark/Benchmark_visitor.hpp>
#include <CGAL/Benchmark/benchmark_format.hpp>

CGAL_BENCHMARK_BEGIN_NAMESPACE

static Benchmark_visitor* visitor; /* global visitor used during parsing */

// Public parser interface, check with decl. in benchmark_parser.h
// ---------------------------------------------------------------

// Opens file 'name' and parses it. Uses visitor 'v' while parsing.
// Returns false if something went wrong. See the visitor for details.
bool benchmark_parse_file( std::string name, Benchmark_visitor* v);

// Starts parsing from stream 'in' with the associated filename 'name'
// (or analogous meaning for different streams) counting linenumbers
// starting from 'n'. Uses visitor 'v' while parsing. Returns false if
// something went wrong. See the visitor for details of the error reporting.
bool benchmark_parse_stream( std::istream&  in, std::string name,
                             Benchmark_visitor* v, int n);

CGAL_BENCHMARK_END_NAMESPACE

using namespace CGAL::benchmark;

/* declaration for flex parser call yylex */
int yylex( void);

/* error function called for parse errors */
void yyerror( const char *s) { visitor->parse_error( std::string(s)); }

static int int_sequence_counter, polynom_varcount;
//static std::string numbertype_name;

struct Number_type {
  virtual ~Number_type() {}
  virtual bool is_equal( Number_type *other ) = 0;
  virtual std::string name() = 0;
};

struct Integer_type : public Number_type {
  virtual bool is_equal( Number_type *other ) {
    return dynamic_cast<Integer_type*>( other ) != NULL;
  }
  virtual std::string name() { return "Integer"; }
};

struct Rational_type : public Number_type {
  virtual bool is_equal( Number_type *other ) {
    return dynamic_cast<Rational_type*>( other ) != NULL;
  }
  virtual std::string name() { return "Rational"; }
};

struct Sqrt_extension_type 
 : public Number_type 
{
  Number_type *t0;
  Number_type *t1;

  Sqrt_extension_type( Number_type *t0, Number_type *t1 ) : t0(t0), t1(t1) {}
  
  virtual ~Sqrt_extension_type() {
    delete t0;
    delete t1;
  }

  virtual bool is_equal( Number_type *other ) {
    Sqrt_extension_type *other_ext = dynamic_cast<Sqrt_extension_type*>( other );
    return t0->is_equal( other_ext->t0 ) &&
           t1->is_equal( other_ext->t1 );
  }

  virtual std::string name() { return "Sqrt_extension<" + t0->name() + "," + t1->name() + ">"; }
};

static std::stack< Number_type* > number_type_stack;
static Number_type *expected_number_type = NULL;


/* Use C++ std::string as semantic value to communicate with lexer */
#define YYSTYPE std::string

%}

/* Elementary data types */
/* --------------------- */
%token INTEGER
%token FNUMBER
%token STRING
%token ERROR
%token UNKNOWN_TOKEN
%token MINUS_INFTY
%token PLUS_INFTY
%token COUNTERCLOCKWISE
%token CLOCKWISE
%token VOID
%token ST GT

/* Structure tokens */
/* ---------------- */
%token FileFormat
%token BenchmarkName
%token Classification

%token List
%token Polynomial_1
%token Polynomial
%token Point_2
%token Degree
%token Vector_2

%token ConicPoint_2

%token LineSegment_2
%token Polyline_2

%token Rotate_2
%token Translate_2

%token Circle_2
%token CircleArc_2

%token IsoEllipse_2
%token IsoEllipseArc_2

%token Conic_2
%token ConicArc_2

%token Cubic_2

%token Quadric_3

%token CSGOperationUnion_3
%token CSGOperationIntersection_3
%token CSGOperationNegation_3

%token Integer
%token Rational
%token Sqrt_extension
%token AlgebraicReal

%%
/* Grammar */
/* ======= */

input: /* an input is a potentially empty sequence of files */
    /* empty */
  | error                            { /* parse error restart here */ }
  | input file
  ;

file: /* must start with unique header */
    file_format file_header_options file_classification file_body
  ;

file_format: /* mandatory fileformat descriptor */
    error_rules                      {}
  | FileFormat '(' STRING ',' INTEGER ',' INTEGER ')'
      { visitor->accept_file_format( $3, atoi( $5.c_str()), atoi( $7.c_str()),
                                     std::string("")); }
  | FileFormat '(' STRING ',' INTEGER ',' INTEGER ',' STRING ')'
      { visitor->accept_file_format( $3, atoi( $5.c_str()), atoi( $7.c_str()), 
                                     $9); }
  ;


file_header_options: /* sequence of optional file header entries */
    /* empty */
  | file_header_options file_header_option
  ;

file_header_option: /* single optional file header entry */
    BenchmarkName '(' STRING ')'     { visitor->accept_benchmark_name( $3); }
  ;

file_classification: /* */
    error_rules             {}
  | classificat
  | file_classification classificat
  ;

classificat: /* */
    Classification '(' STRING ',' STRING ',' STRING ',' 
                       STRING ',' STRING ',' STRING ')'
      { visitor->accept_classification( $3, $5, $7, $9, $11, $13); }
  ;

file_body: /* sequence of statements (stmt) */
    /* empty */
  | file_body stmt
  ;

stmt_sequence: /* comma separated sequence of statements (stmt) */
    /* empty */
  | stmt_sequence_non_empty
  ;

stmt_sequence_non_empty: /* comma separated sequence of statements (stmt) */
    stmt
  | stmt_sequence_non_empty ',' stmt
  ;

stmt: /* */
    error                            { /* parse error restart here */ }
    error_rules                      {}
  | List                             { visitor->begin_list(); }
      '(' stmt_sequence ')'          { visitor->end_list(); }
  | LineSegment_2 line_segment_2
  | Polyline_2 polyline_2
  | Circle_2 circle_2
  | CircleArc_2 circle_arc_2
  | IsoEllipse_2 iso_ellipse_2
  | IsoEllipseArc_2 iso_ellipse_arc_2
  | Conic_2 conic_2
  | ConicArc_2 conic_arc_2
  | transformed_conic
  | Cubic_2 '(' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ','
                INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ')'
      { visitor->accept_cubic_2($3, $5, $7, $9, $11, $13, $15,$17, $19, $21); }
  | Quadric_3 quadric_3
  | csg_operation_3
  | polynomial
  ;

line_segment_2:
    '('                              { visitor->begin_line_segment_2(); }
      line_segment_2_body ')'        { visitor->end_line_segment_2(); }
  ;

line_segment_2_body:
    point_2 ',' point_2
  ;

polyline_2:
    '('                              { visitor->begin_polyline_2(); }
      polyline_2_body ')'            { visitor->end_polyline_2(); }
  ;

polyline_2_body:
    point_2 ',' point_2_list
  ;

circle_2:
    '('                              { visitor->begin_circle_2(); }
      circle_2_body ')'              { visitor->end_circle_2(); }
  ;

circle_2_body:
    point_2 ',' intorrat
  | point_2 ',' point_2 ',' point_2
  ;

circle_arc_2:
    '('                              { visitor->begin_circle_arc_2(); }
      circle_arc_2_body ')'          { visitor->end_circle_arc_2(); }    

circle_arc_2_body:
    point_2 ',' point_2 ',' point_2
  ;

iso_ellipse_2:
    '('                              { visitor->begin_iso_ellipse_2(); }
      iso_ellipse_2_body ')'         { visitor->end_iso_ellipse_2(); }
  ;

iso_ellipse_2_body:
    point_2 ',' intorrat ',' intorrat
  ;

iso_ellipse_arc_2:
    '('                              { visitor->begin_iso_ellipse_arc_2(); }
      iso_ellipse_arc_2_body ')'     { visitor->end_iso_ellipse_arc_2(); }
  ;

iso_ellipse_arc_2_body:
    IsoEllipse_2 iso_ellipse_2 ',' point_2 ',' point_2 ','
      orientation                    { visitor->accept_orientation( $8); }
  ;

polynomial: /* */
    // directly encoding '<' does not work, interestingly ...
    Polynomial ST INTEGER optional_typename GT
      { 
        polynom_varcount = atoi( $3.c_str() );
        visitor->begin_polynomial( polynom_varcount,
        expected_number_type->name() );
      }
      '(' monom_list ')'             {
                                       visitor->end_polynomial(); 
                                       delete expected_number_type; 
                                     }
  | Polynomial_1                     { visitor->begin_polynomial_1(); }
      '(' integer_sequence1 ')'      { visitor->end_polynomial_1(); }
  ;

optional_typename :
    /* empty */              { expected_number_type = new Integer_type(); }
  | ',' typename             { 
                               expected_number_type = number_type_stack.top();
                               number_type_stack.pop(); 
                               assert( number_type_stack.size() == 0 );
                             }
  ;

typename :
    Sqrt_extension ST typename ',' typename GT 
      { 
        assert( number_type_stack.size() >= 2 );
        Number_type *t1 = number_type_stack.top(); 
        number_type_stack.pop();
        Number_type *t0 = number_type_stack.top();
        number_type_stack.pop();
        number_type_stack.push( new Sqrt_extension_type( t0, t1 ) ); 
      }
  |  Rational                { number_type_stack.push( new Rational_type() ); }
  |  Integer                 { number_type_stack.push( new Integer_type() ); }
  ;  

monom_list:
    monom
  | monom_list ',' monom
  ;

point_2_list:
    point_2
  | point_2_list ',' point_2
  ;

monom:
    '(' numeric_value ','
      {
        visitor->begin_monom( $2 );
        int_sequence_counter = 0;
      }
    '(' integer_sequence1 ')' ')'
      {
        visitor->end_monom();
        if( int_sequence_counter != polynom_varcount ) {
          std::stringstream errmsg;
          errmsg << "expected monom with " << polynom_varcount 
		 << " variables.";
          yyerror( errmsg.str().c_str() );
        }
      }
  ;

numeric_value:
    numeric_value1
      { 
        Number_type *actual_number_type = number_type_stack.top();
        number_type_stack.pop();
        if( ! expected_number_type->is_equal( actual_number_type ) ) {
          yyerror( ("type mismatch. expected \n" + 
                    expected_number_type->name() + 
                    "\nbut got\n" + actual_number_type->name()).c_str() );
          exit(EXIT_FAILURE);
        } else
          $$ = $1; 
        delete actual_number_type;
      }
  ;
  
  
numeric_value1:
    INTEGER
      {
        $$ = $1;
        number_type_stack.push( new Integer_type() );
      }
  | Rational '(' INTEGER ',' INTEGER ')' 
      {
        $$ = "RAT[" + $3 + "," + $5 + "]";
        number_type_stack.push( new Rational_type() );
      }
  | Sqrt_extension '(' numeric_value1 ',' numeric_value1 ',' numeric_value1 ')'
      {
        $$ = "EXT[" + $3 + "," + $5 + "," + $7 + "]";
        assert( number_type_stack.size() >= 3 );
        Number_type *t2 = number_type_stack.top(); number_type_stack.pop();
        Number_type *t1 = number_type_stack.top(); number_type_stack.pop();
        Number_type *t0 = number_type_stack.top(); number_type_stack.pop();

        if( ! t1->is_equal( t0 ) ) {
          yyerror( "type mismatch. for sqrt_ext, t0 and t1 should be equal" );
          exit(EXIT_FAILURE);
        }
        number_type_stack.push( new Sqrt_extension_type( t0, t2 ) );
      }
  ;

conic: /* */
    error_rules                      {}
  | Conic_2 conic_2
  | transformed_conic
  ;

conic_2: /* */
    error_rules                      {}
  | '(' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ')'
      { visitor->accept_conic_2($2, $4, $6, $8, $10, $12); }
  ;


transformed_conic:  /* */
    Rotate_2 '('                     { visitor->begin_rotate_2(); }
      Degree '(' INTEGER ')'         { visitor->accept_integer( $6); }
      ',' point_2 ',' conic ')'      { visitor->end_rotate_2(); }
  | Translate_2 '('                  { visitor->begin_translate_2(); }
      Vector_2 '('                   { visitor->begin_vector_2(); }
        rational ',' rational ')'    { visitor->end_vector_2(); }
      ',' conic ')'                  { visitor->end_translate_2(); }
  ;

conic_arc_2:
    '('                              { visitor->begin_conic_arc_2(); }
      conic_arc_2_body ')'           { visitor->end_conic_arc_2(); }
  ;

conic_arc_2_body: /* */
     error_rules {}
  | Conic_2 conic_2 ',' 
      ConicPoint_2 conic_point_2 ',' ConicPoint_2 conic_point_2 ','
      orientation                    { visitor->accept_orientation( $10); }
  | Conic_2 conic_2 ',' point_2 ',' point_2 ','
      orientation                    { visitor->accept_orientation( $8); }

  | ConicPoint_2 conic_point_2
  ;

integer_sequence1: /* comma separated integers */
    INTEGER          
      { visitor->accept_integer( $1); ++int_sequence_counter; }
  | integer_sequence1 ',' INTEGER 
      { visitor->accept_integer( $3); ++int_sequence_counter; }
  ;

point_2: /* a 2d point */
    Point_2 '(' INTEGER ',' INTEGER ')' 
                                     { visitor->accept_point_2( $3, $5); }
  | Point_2 '(' INTEGER ',' INTEGER ',' INTEGER ')'
                                     { visitor->accept_point_2( $3, $5, $7); }
  | Point_2 '(' Rational '(' INTEGER ',' INTEGER ')' ','
                Rational '(' INTEGER ',' INTEGER ')' ')'
                                { visitor->accept_point_2( $5, $7, $12, $14); }
  ;

conic_point_2: /* */
    '('                              { visitor->begin_conic_point_2();}
    Conic_2 conic_2 ',' algorint ')' { visitor->end_conic_point_2(); }
  ;

algorint: /* */
    algebraic_real ',' inti
  | infty                            { visitor->accept_infty( $1); }
      ',' INTEGER                    { visitor->accept_integer( $3); }
  ;


algebraic_real: /* */
    AlgebraicReal '('                { visitor->begin_algebraic_real(); }
      Polynomial_1 '('               { visitor->begin_polynomial_1(); }
        integer_sequence1 ')'        { visitor->end_polynomial_1(); }
      ',' rational ',' rational ',' 
      INTEGER                        { visitor->accept_integer($15); }
    ')'                              { visitor->end_algebraic_real(); }
  ;

intorrat: /* */
    INTEGER                          { visitor->accept_integer( $1); }
  | rational
  ;

rational: /* */
   Rational '(' INTEGER ',' INTEGER ')'
                                     { visitor->accept_rational( $3, $5); }
  ;

inti: /* */
    algebraic_real
  | INTEGER                          { visitor->accept_integer( $1); }
  | infty                            { visitor->accept_infty( $1); }
  ;

infty: /* */
    MINUS_INFTY
  | PLUS_INFTY
  ;

orientation: /* */
    COUNTERCLOCKWISE
  | CLOCKWISE
  | VOID
  ;

/*
double_val: // a double value can be either an FNUMBER or an INTEGER
    FNUMBER
  | INTEGER                          { visitor->accept_integer( $1); }
  ;
*/

error_rules: /* */
    ERROR                            { visitor->parse_error( $1); }
  | UNKNOWN_TOKEN                    { visitor->unknown_token( $1); }
  ;

quadric_3:
    '(' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ','
        INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ')'
      { visitor->accept_quadric_3( $3, $5, $7, $9, $11, 
                                   $13, $15, $17, $19, $21); }
  ;

csg_operation_3:
    CSGOperationUnion_3              { visitor->begin_csg_union_3(); }
      csg_operation_union_3          { visitor->end_csg_union_3(); }
  | CSGOperationIntersection_3       { visitor->begin_csg_intersection_3(); }
      csg_operation_intersection_3   { visitor->end_csg_intersection_3(); }
  | CSGOperationNegation_3           { visitor->begin_csg_negation_3(); }
      csg_operation_negation_3       { visitor->end_csg_negation_3(); }
  ;

csg_operation_union_3:
    '(' quadric_3 ',' quadric_3 ')'
  | '(' quadric_3 ',' csg_operation_3 ')'
  | '(' csg_operation_3 ',' quadric_3 ')'
  | '(' csg_operation_3 ',' csg_operation_3 ')'
  ;

csg_operation_intersection_3:
    '(' quadric_3 ',' quadric_3 ')'
  | '(' quadric_3 ',' csg_operation_3 ')'
  | '(' csg_operation_3 ',' quadric_3 ')'
  | '(' csg_operation_3 ',' csg_operation_3 ')'
  ;

csg_operation_negation_3:
    '(' quadric_3 ')'
  | '(' csg_operation_3 ')'
  ;

/* End of Grammar */
/* ============== */
%%

CGAL_BENCHMARK_BEGIN_NAMESPACE

// Opens file 'name' and parses it. Uses visitor 'v' while parsing.
// Returns false if something went wrong. See the visitor for details.
bool benchmark_parse_file( std::string name, Benchmark_visitor* v) {
    std::ifstream in( name.c_str());
    if ( ! in) {
        v->parse_error( std::string( "cannot open file '") + name +
                        std::string( "'."));
        return false;
    }
    visitor = v;
    benchmark_init_lexer( in, name, 1);
    yyparse();
    return ! v->error();
}

// Starts parsing from stream 'in' with the associated filename 'name'
// (or analogous meaning for different streams) counting linenumbers
// starting from 'n'. Uses visitor 'v' while parsing. Returns false if
// something went wrong. See the visitor for details of the error reporting.
bool benchmark_parse_stream( std::istream&  in, std::string name,
                             Benchmark_visitor* v, int n) {
    visitor = v;
    benchmark_init_lexer( in, name, n);
    yyparse();
    return ! v->error();
}

CGAL_BENCHMARK_END_NAMESPACE

/* EOF */
