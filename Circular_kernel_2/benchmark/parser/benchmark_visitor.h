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
**************************************************************************/

#ifndef BENCHMARK_VISITOR_H
#define BENCHMARK_VISITOR_H

#include <string>

// Exception thrown by Benchmark_visitor if mode is THROW_ERRORS
struct Benchmark_parser_exception {
    std::string message;
    Benchmark_parser_exception() {}
    Benchmark_parser_exception( std::string s) : message(s) {}
};

// Visitor pattern for the benchmark parser
class Benchmark_visitor {
public:
    enum Error_mode { IGNORE_ERRORS = 0, PRINT_ERRORS = 1, 
                      THROW_ERRORS  = 2, EXIT_ERRORS  = 4 };
protected:
    // error handling
    int          m_mode;          // Error_mode or'ed together
    bool         m_error;         // flags if an error has occured
    std::string  m_error_message; // error message kept from last error
    virtual void error_handler( std::string s);

    // information from the benchmark file already understood and parsed
    int          m_major;
    int          m_minor;
    std::string  m_format_options;
    std::string  m_benchmark_name;

public:    
    Benchmark_visitor( int mode = PRINT_ERRORS )
        : m_mode(mode), m_error( false), m_error_message(""), 
          m_major(-1), m_minor(-1) {}
    virtual ~Benchmark_visitor() {}

    // error status and mode
    bool         error() const { return m_error; }
    std::string  error_message() const { return m_error_message; }
    void         reset_error() {
        m_error = false;
        m_error_message = std::string( "");
    }

    // information from the benchmark file already understood and parsed
    int          version_major()  const { return m_major; }
    int          version_minor()  const { return m_minor; }

    virtual void reset_header();

    // Options are comma separated. We have also a comma left and one right 
    // to have simple contains-an-option ",OPT," queries.
    std::string  format_options() const { return m_format_options; }
    std::string  benchmark_name() const { return m_benchmark_name; }

    // error handlers
    virtual void parse_error( std::string s);
    virtual void unknown_token( std::string s);
    virtual void token_not_handled( std::string s);
    void tnh( std::string s) { token_not_handled(s); } // shortcut

    // file header entries
    virtual void accept_file_format( std::string s, int major, int minor,
                                     std::string options);
    virtual void accept_benchmark_name( std::string s);
    virtual void accept_classification( std::string problem,
                                        std::string geom, 
                                        std::string clas,
                                        std::string family,
                                        std::string instance,
                                        std::string release);

    // terminal nodes
    virtual void accept_infty( std::string s) { tnh( "Infty"); }
    virtual void accept_orientation( std::string s) { tnh( "Orientation"); }

    virtual void accept_integer( std::string s) { tnh( "Integer"); }
    virtual void accept_rational( std::string num, std::string denom) { 
        tnh( "Rational"); 
    }
    virtual void accept_fnumber( double d) { tnh( "Fnumber"); }
    virtual void accept_string( std::string s) { tnh( "String"); }
    virtual void accept_minus_infty (std::string s) { tnh("Minus_infty"); }
    virtual void accept_plus_infty (std::string s) { tnh("Plus_infty"); }
    virtual void accept_counter (std::string s) { tnh("Counterclockwise"); }
    virtual void accept_clockwise (std::string s) { tnh("Clockwise"); }
    virtual void accept_void (std::string s) { tnh("Void"); }

    virtual void accept_point_2( std::string x, std::string y) { 
        tnh( "Point_2(x,y)"); 
    }
    virtual void accept_point_2( std::string x, std::string y, std::string w) {
        tnh( "Point_2(x,y,w)"); 
    }
    virtual void accept_point_2( std::string x_num, std::string x_denom, 
                                 std::string y_num, std::string y_denom) {
        tnh( "Point_2(x_num, x_denom, y_num, y_denom)"); 
    }

    virtual void begin_algebraic_real() { tnh( "Begin_algebraic_real" ); }
    virtual void end_algebraic_real()   { tnh( "End_algebraic_real" ); }


    virtual void begin_conic_point_2() {
        tnh( "Begin_conic_point_2" );
    }

    virtual void end_conic_point_2() {
        tnh( "End_conic_point_2" );
    }
    
    // non-terminal nodes
    virtual void begin_list() { tnh( "Begin_list"); }
    virtual void end_list() { tnh( "End_list"); }

    virtual void begin_polynomial_1() { tnh( "Begin_polynomial_1"); }
    virtual void end_polynomial_1()   { tnh( "End_polynomial_1"); }

    virtual void begin_line_segment_2() { tnh( "Begin_line_segment_2"); }
    virtual void end_line_segment_2() { tnh( "End_line_segment_2"); }

    virtual void accept_conic_2(std::string A, std::string B, std::string C,
                                std::string D, std::string E, std::string F) { 
        tnh("Conic_2"); 
    }


    virtual void begin_conic_arc_2() { tnh( "Begin_conic_arc_2" );}
    virtual void end_conic_arc_2()   { tnh( "End_conic_arc_2" );}
    
    virtual void begin_circle_2() { tnh("Begin_circle_2"); }
    virtual void end_circle_2() { tnh("End_circle_2"); }
 
    virtual void accept_cubic_2(std::string A, std::string B, std::string C,
                                std::string D, std::string E, std::string F,
                                std::string G, std::string H, std::string K,
                                std::string L) {
        tnh("Cubic_2");
    }


    virtual void accept_quadric_3(std::string A, std::string B, std::string C,
                                  std::string D, std::string E, std::string F,
                                  std::string G, std::string H, std::string K,
                                  std::string L) {
        tnh("Quadric_3");
    }
 virtual void begin_CircularArc_2(){
        tnh("begin_CircularArc_2");
        }
    virtual void end_CircularArc_2(){
        tnh("end_CircularArc_2");
        }
   virtual void begin_LineArc_2(){
        tnh("begin_LineArc_2");
        }
    virtual void end_LineArc_2(){
        tnh("end_LineArc_2");
        }
   virtual void begin_CircularPoint_2(){
        tnh("begin_CircularPoint_2");
        }
    virtual void end_CircularPoint_2(){
        tnh("end_CircularPoint_2");
        }
};


#endif // BENCHMARK_VISITOR_H //
