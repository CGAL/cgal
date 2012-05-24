/**************************************************************************
// Copyright (c) 2004  Max-Planck-Institut Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of BenchmarkParser; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s) : Lutz Kettner
//             Franziska Ebert <febert@mpi-sb.mpg.de>
**************************************************************************/

%{
/* C/C++ declaration section */
/* ========================= */
//#include <cstdlib>  /* for atoi */
#include <string>    /* for std::string */
#include <fstream>   /* for std::ifstream */
#include <benchmark_parser.h>
#include <benchmark_visitor.h>
#include <iostream>
static Benchmark_visitor* visitor; /* global visitor used during parsing */

/* declaration for flex parser call yylex */
int yylex( void);

/* error function called for parse errors */
void yyerror( char *s) { visitor->parse_error( std::string(s)); }

/* Use C++ std::string as semantic value to communicate with lexer */
#define YYSTYPE std::string

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

/* Structure tokens */
/* ---------------- */
%token FileFormat
%token BenchmarkName
%token Classification

%token List
%token Rational
%token Polynomial_1
%token Point_2
%token AlgebraicReal


%token ConicPoint_2
%token LineSegment_2
%token Conic_2

%token CircularArc_2
%token LineArc_2
%token CircularPoint_2

%token ConicArc_2
%token Circle_2
%token Cubic_2


%token Quadric_3


%%
/* Grammar */
/* ======= */

input: /* an input is a potentially empty sequence of files */
  /* empty */
  | error                       { /* parse error restart here */ }
  | input file
;

file: /* must start with unique header */
    file_format file_header_options file_classification file_body
;

file_format: /* mandatory fileformat descriptor */
    error_rules                 {}
  | FileFormat '(' STRING ',' INTEGER ',' INTEGER ')' 
                                { visitor->accept_file_format( $3, 
                                      atoi( $5.c_str()), atoi( $7.c_str()),
                                      std::string("")); }
  | FileFormat '(' STRING ',' INTEGER ',' INTEGER ',' STRING ')' 
                                { visitor->accept_file_format( $3, 
                                     atoi( $5.c_str()), atoi( $7.c_str()),
                                                              $9); }    
;

file_header_options: /* sequence of optional file header entries */
    /* empty */
  | file_header_options file_header_option
;

file_header_option: /* single optional file header entry */
    BenchmarkName '(' STRING ')' { visitor->accept_benchmark_name( $3); } 
;

file_classification: /* */
    error_rules             {}
  | classificat 
  | file_classification classificat 
;

classificat: /* */
  Classification '(' STRING ',' STRING ',' STRING ',' STRING ',' 
                     STRING ',' STRING ')'
                            { visitor->accept_classification( $3, $5, $7, 
                                                              $9, $11, $13); }
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

stmt:
    error                        { /* parse error restart here */ }
    error_rules                  {}
  | List                         { visitor->begin_list(); }
       '(' stmt_sequence ')'     { visitor->end_list(); }
  | Circle_2      		 { visitor->begin_circle_2(); } 
       '(' circle_2 ')'   	 { visitor->end_circle_2(); }
  | Conic_2 conic_2
  | LineSegment_2                { visitor->begin_line_segment_2(); } 
       '(' point_2 ',' point_2 ')' 
                                 { visitor->end_line_segment_2(); } 
  | ConicArc_2 conic_arc_2 

  | Cubic_2 '(' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' 
                INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ')'
                                 { visitor->accept_cubic_2( $3, $5, $7, $9, 
                                                            $11, $13, $15, 
                                                            $17, $19, $21); }
  | Quadric_3 '(' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' 
                  INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ',' INTEGER ')'
                                 { visitor->accept_quadric_3( $3, $5, $7, $9, 
                                                              $11, $13, $15, 
                                                              $17, $19, $21); }

  | CircularPoint_2 		{visitor->begin_CircularPoint_2();}
	'('circular_arc_point')'
				{visitor->end_CircularPoint_2();}
  | LineArc_2			{visitor->begin_LineArc_2();}
	'('line_arc_2')'	
				{visitor->end_LineArc_2();}
  |CircularArc_2         	{ visitor->begin_CircularArc_2();} 
      '('circular_arc_2 ')'  
                                 { visitor->end_CircularArc_2();} 
;
circular_arc_point:
circular_arc_point:
error_rules{}
 | point_2
 | AlgebraicReal ',' AlgebraicReal
;
line_arc_2:
error_rules {}
 | LineSegment_2
 | point_2 ',' point_2 
 | CircularPoint_2 '(' circular_arc_point ')' 
       ',' CircularPoint_2 '(' circular_arc_point ')'
;
circular_arc_2:
error_rules {}
 | Circle_2 '(' circle_2 ')'
 | Circle_2 '(' circle_2 ')' ',' CircularPoint_2 '(' circular_arc_point ')' 
       ',' CircularPoint_2 '(' circular_arc_point ')'
 | point_2 ',' point_2 ',' rational
;
circle_2:
error_rules {}
|  point_2 ',' rational  
|  point_2 ',' INTEGER  { visitor->accept_integer($3); } 

;
conic_2: /* */
    error_rules                  {}

  | '(' INTEGER ',' INTEGER ',' INTEGER ',' 
        INTEGER ',' INTEGER ',' INTEGER ')' 
                                 { visitor->accept_conic_2( $2, $4, $6, 
                                                            $8, $10, $12); }


;
conic_arc_2: /* */
     error_rules {}
  | '(' Conic_2                  { visitor->begin_conic_arc_2(); } 
            conic_2
        ',' ConicPoint_2 conic_point_2 
        ',' ConicPoint_2 conic_point_2
	',' orientation          { visitor->accept_orientation( $12); }
     ')' 
                                 { visitor->end_conic_arc_2(); }


  | '(' ConicPoint_2             { visitor->begin_conic_arc_2(); }
            conic_point_2 ')'
                                 { visitor->end_conic_arc_2(); }
            
;

integer_sequence1: /* comma separated integers */
    error_rules                   {}
  | INTEGER                       { visitor->accept_integer( $1); }
  | integer_sequence1 ',' INTEGER { visitor->accept_integer( $3); }
;


point_2: /* a 2d point */
    error_rules                   {}
  | Point_2 '(' INTEGER ',' INTEGER ')' 
                                  { visitor->accept_point_2( $3, $5); }
  | Point_2 '(' INTEGER ',' INTEGER ',' INTEGER ')' 
                                  { visitor->accept_point_2( $3, $5, $7); }
  | Point_2 '(' Rational '(' INTEGER ',' INTEGER ')' ','
                Rational '(' INTEGER ',' INTEGER ')' ')' 
                                  { visitor->accept_point_2( $5, $7, 
                                                             $12, $14); }
;

conic_point_2: /* */
     error_rules                  {}
  | '(' Conic_2                   { visitor->begin_conic_point_2();} 
            conic_2
        ',' algorint
     ')'                          { visitor->end_conic_point_2(); }

;

algorint: /* */
    error_rules                   {}
  | algebraic_real ',' inti
  | infty ',' INTEGER             { visitor->accept_infty( $1); 
                                    visitor->accept_integer( $3); }
;


algebraic_real: /* */
    error_rules                   {}
  | AlgebraicReal                 { visitor->begin_algebraic_real(); }
        '(' Polynomial_1          { visitor->begin_polynomial_1(); }
              '(' integer_sequence1 ')' 
                                  { visitor->end_polynomial_1(); }

        ',' rational
        ',' rational
        ',' INTEGER               { visitor->accept_integer($15); }
        ')'                       { visitor->end_algebraic_real(); }
;

rational: /* */
    error_rules                   {}
  | Rational '(' INTEGER ',' INTEGER ')' 
                                  { visitor->accept_rational( $3, $5); }
;

inti: /* */ 
    INTEGER                       { visitor->accept_integer( $1); }  
  | infty                         { visitor->accept_infty( $1); }
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

double_val: /* a double value can be either an FNUMBER or an INTEGER */
    FNUMBER
  | INTEGER                     { visitor->accept_integer( $1); } 
;

error_rules:
    ERROR                       { visitor->parse_error( $1); }
  | UNKNOWN_TOKEN               { visitor->unknown_token( $1); }
;

/* End of Grammar */
/* ============== */
%%

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

/* EOF */
