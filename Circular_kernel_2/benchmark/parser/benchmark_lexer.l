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

%{
#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>

/* Use C++ std::string as semantic value to communicate with lexer */
#define YYSTYPE std::string
#include <benchmark_parser.tab.h> // gets definitions from syntax parser

// Public lexer interface, check with decl. in benchmark_parser.h
// --------------------------------------------------------------

// Current input stream and associated values.
std::istream* benchmark_in         = & std::cin;
int           benchmark_linenumber = 1;
std::string   benchmark_filename(    "<cin>");

// Initialize lexer to scan from input stream with name and linenumber.
// The caller is responsible for the lifetime of 'in' that must live
// during the scan.
void benchmark_init_lexer( 
    std::istream& in, std::string name, int linenumber = 1);

// Reset lexer to clean state.
void benchmark_reset_lexer();

// Writes a trace of the current include file nesting to the 'out' stream.
// Appends the 'fill' string after each file listed.
void benchmark_include_file_trace( std::ostream& out, std::string fill);


// Private lexer interface, all static to keep it local
// ----------------------------------------------------

// Stack of input streams: the item stored in the stack.
struct Include_stack_item {
    std::istream*    in;         // input stream
    int              linenumber; // linenumber in file for error messages
    std::string      filename;   // filename (or similar string)
    YY_BUFFER_STATE  buffer;     // Opaque pointer to flex-buffer struct.
    Include_stack_item() {}
    Include_stack_item( std::istream* i, int ln, std::string name,
                        YY_BUFFER_STATE  buf)
        : in(i), linenumber(ln), filename(name), buffer(buf) {}
};

// Stack of input streams is a std::list
typedef std::list<Include_stack_item> Include_stack;

// We have one static global variable for the stack and some static functions.
// The include stack is except for the trace() debug function not 
// accessible from other files (is not polluting namespaces) since 
// the Include() function and the <<EOF>> condition are all handled here.
// The function implementations are at the end of this file.
static Include_stack include_stack;

// Push current state on include_stack, initializes lexer with new file
// input stream, name, and linenumber.
static void open_include_stream( 
    std::istream *in, std::string name, int linenumber = 1);

// Shortcut for open file with name and then call open_include_stream.
// Returns false if the file could not be opened successfully.
static bool open_include_file( std::string name, int linenumber = 1);

// Closes current input stream. Initializes lexer with new file input stream,
// name, and linenumber, from the include_stack's top. Pop's an element from
// the include_stack. Returns false if include_stack was empty.
static bool close_include_stream();


// Make the lexer read from benchmark_in instead of yyin
#define YY_INPUT(buf,result,max_size) { \
    benchmark_in->read( buf, max_size);  \
    result = benchmark_in->gcount();     \
}

// Count newlines in the string s
static void count_newlines( const char* s) {
    while ( *s) {
        if ( *s == '\n')
            ++benchmark_linenumber;
        ++s;
    }
}

static int  comment_nesting = 0; // counts nesting depth of () in Comments

#define YY_BREAK  /* a do nothing */

/* --------------------------------------------------------------------
    Parsing Modes:
      --  INITIAL:      main mode for sequence of tokens
      --  IncludeMode:  parses lciInclude filename,
      --  CommentMode:  Comment(...) parsing of nested parantheses
                        # comments and strings are correctly ignored
-------------------------------------------------------------------- */
%}

%x IncludeMode
%x CommentMode

space           [\t \r]*
newline         [\n]
spcnl           [\t \n\r]*
integer         [+-]?[0-9]+
/* An fnumber consists of up to three parts:
 *  1) a sign [+-]
 *  2) a mantissa (digits [0-9] with at most one dot in between) 
 *  3) an exponent [eE]{integer}
 *
 *  Part 1) is optional.
 *  Part 2) must contain at least one digit,
 *          but may have leading or trailing dot
 *  Part 3) is optional iff the mantissa contains a dot
 *          (so that fnumbers and integers remain distinguishable).
 */
fnumber [+-]?(([0-9]+[.][0-9]*|[.][0-9]+)([eE]{integer})?|[0-9]+[eE]{integer})
idfier          [a-zA-Z][a-zA-Z0-9_]*
minus_infty     MINUS_INFTY
plus_infty      PLUS_INFTY
filename        [^ \t\n\r\\\{\}\[\]()]+
void            VOID
counterclockwise COUNTERCLOCKWISE
clockwise       CLOCKWISE

%%
 /* Tokens    */
 /* --------- */
"FileFormat"                    { return FileFormat;}
"BenchmarkName"                 { return BenchmarkName; }
"Classification"                { return Classification; }
"List"                          { return List; }
"Rational"                      { return Rational; }
"Polynomial_1"                  { return Polynomial_1; }
"Point_2"                       { return Point_2; }
"AlgebraicReal"                 { return AlgebraicReal; }


"ConicPoint_2"                  { return ConicPoint_2; } 
"Conic_2"                       { return Conic_2; }


"ConicArc_2"                    { return ConicArc_2; }
"Circle_2"                      { return Circle_2; }
"LineSegment_2"                 { return LineSegment_2; }
"Cubic_2"                       { return Cubic_2; }
"CircularArc_2"			{ return CircularArc_2;}
"LineArc_2"			{ return LineArc_2;}
"CircularPoint_2"		{return CircularPoint_2;}
"Quadric_3"                     { return Quadric_3; }


 /* Integer numbers and float numbers */
 /* --------------------------------- */
{integer}                       { yylval = std::string( yytext); 
                                  return INTEGER; 
                                }
{fnumber}                       { yylval = std::string( yytext); 
                                  return FNUMBER; 
                                }
{minus_infty}                   { yylval = std::string( yytext);
                                  return MINUS_INFTY;
                                }
{plus_infty}                    { yylval = std::string( yytext);
                                  return PLUS_INFTY;
                                }
{void}                          { yylval = std::string( yytext);
                                  return VOID;
                                }
{counterclockwise}              { yylval = std::string( yytext);
                                  return COUNTERCLOCKWISE;
                                }
{clockwise}                     { yylval = std::string( yytext);
                                  return CLOCKWISE;
                                }


 /* Non-recognized tokens are errors */
 /* -------------------------------- */
{idfier}                        { yylval = std::string( yytext); 
                                  return UNKNOWN_TOKEN; 
                                }

 /* Handle include files      */
 /* ------------------------- */
"IncludeFile"{spcnl}[(]{spcnl}  { BEGIN( IncludeMode);
                                  count_newlines( yytext);
                                  break;
                                }
<IncludeMode>{filename}         { /* remove remaining chars before the ')' */
                                  int c = yyinput();
                                  while( c != EOF && c != ')') {
                                      if ( c == '\n')
                                          ++benchmark_linenumber;
                                      c = yyinput();
                                  }
                                  BEGIN( INITIAL);
                                  if ( c == EOF) {
                                      yylval = std::string( "Reached EOF while"
                                          " parsing include filename '") + 
                                          std::string( yytext) + 
                                          std::string( "'.");
                                      return ERROR;
                                  }
                                  open_include_file( yytext);
                                  break;
                                }

 /* Strings with quoted \" and \\  */
 /* ------------------------------ */
<INITIAL,CommentMode>"\""       { int c = yyinput();
                                  yylval = std::string("");
                                  bool quoted_char = false; 
                                  while (c != EOF && (c!='"' || quoted_char)){
                                      if ( c == '\n')
                                          ++benchmark_linenumber;
                                      if ( ! quoted_char && c == '\\') {
                                          quoted_char = true;
                                      } else {
                                          quoted_char = false;
                                          yylval.push_back( char(c));
                                      }
                                      c = yyinput();
                                  }
                                  if ( c == EOF) {
                                      yylval = std::string( "Reached EOF while"
                                          " parsing string constant '") + 
                                          yylval + std::string( "'.");
                                      BEGIN( INITIAL);
                                      return ERROR;
                                  }
                                  if ( YY_START == INITIAL) {
                                      return STRING;
                                  }
                                  yylval = std::string("");
                                  break;
                                }

 /* One line comment starting with # */
 /* -------------------------------- */
<INITIAL,CommentMode>"#"        { /* remove remaining chars before the '\n' */
                                  int c = yyinput();
                                  while( c != EOF && c != '\n') {
                                      c = yyinput();
                                  }
                                  if ( c == '\n')
                                      ++benchmark_linenumber;
                                  break;
                                }

 /* Nestable and multi-line Comment( ... ) */
 /* -------------------------------------- */
"Comment"{spcnl}"("{spcnl}      { BEGIN( CommentMode);
                                  comment_nesting = 1;
                                  count_newlines( yytext);
                                  break;
                                }

<CommentMode>{newline}          { ++benchmark_linenumber; break; }

<CommentMode>"("                { ++comment_nesting; break; }

<CommentMode>")"                { if ( --comment_nesting == 0)
                                      BEGIN( INITIAL);
                                  break;
                                }

<CommentMode><<EOF>>            { yylval = std::string( "Reached EOF while "
                                    "parsing 'Comment(...)'.");
                                  BEGIN( INITIAL);
                                  return ERROR;
                                }

<CommentMode>.                  { break; }

 /* Count line numbers for good error messages */
{newline}                       { ++benchmark_linenumber; break; }

 /* Ignore white spaces */
{space}                         { break; }

 /* stop scanning at EOF, maybe continue with surrounding file */
<<EOF>>                         { yyterminate(); }

 /* single characters passed to the parser: */
[(),]                           { return yytext[0]; }

 /* all other single characters are errors */
.                               { yylval = std::string( "Found illegal char '")
                                    + std::string( yytext) + std::string("'.");
                                  return ERROR;
                                }

%%

// implements the proper EOF behavior, pop's file from include stack
extern "C" int yywrap() {
    if ( close_include_stream())
        return 0; // there was a file waiting on the stack, we continue scan
    return 1; // no file in the stack left, we stop scan
}

// Push current state on include_stack, initializes lexer with new file
// input stream, name, and linenumber.
static void open_include_stream( std::istream *in, std::string name, 
                                 int linenumber) 
{
    include_stack.push_front( Include_stack_item( 
        benchmark_in, benchmark_linenumber, benchmark_filename, 
        YY_CURRENT_BUFFER));
    benchmark_in = in;
    benchmark_linenumber = linenumber;
    benchmark_filename = name;
    yy_switch_to_buffer( yy_create_buffer( 0, YY_BUF_SIZE));
}

// Shortcut for open file with name and then call open_include_stream.
// Returns false if the file could not be opened successfully.
static bool open_include_file( std::string name, int linenumber) {
    std::ifstream* in = new std::ifstream( name.c_str());
    if ( ! *in) {
        delete in;
        return false;
    }
    open_include_stream( in, name, linenumber);
    return true;
}

// Closes current input stream. Initializes lexer with new file input stream,
// name, and linenumber, from the include_stack's top. Pop's an element from
// the include_stack. Returns false if include_stack was empty.
static bool close_include_stream() {
    if ( include_stack.empty())
        return false;
    delete benchmark_in; // closes the stream
    benchmark_in         = include_stack.front().in;
    benchmark_linenumber = include_stack.front().linenumber;
    benchmark_filename   = include_stack.front().filename;
    yy_delete_buffer( YY_CURRENT_BUFFER);
    yy_switch_to_buffer( include_stack.front().buffer);
    include_stack.pop_front();
    return true;
}

// Public interface of the lexer component
// ---------------------------------------

// Initialize lexer to scan from input stream with name and linenumber.
// The caller is responsible for the lifetime of 'in' that must last 
// during the scan.
void benchmark_init_lexer( 
    std::istream& in, std::string name, int linenumber) 
{
    benchmark_reset_lexer();
    benchmark_in = & in;
    benchmark_linenumber = linenumber;
    benchmark_filename = name;
    yy_switch_to_buffer( yy_create_buffer( 0, YY_BUF_SIZE));
}

// reset lexer to clean state. Can be used to parse another file then.
void benchmark_reset_lexer() {
    while( close_include_stream()) // close all pre-existing include files
        ; // empty while body
    benchmark_in         = & std::cin;
    benchmark_linenumber = 1;
    benchmark_filename   = std::string("<cin>");
    yyrestart(0);
    BEGIN( INITIAL);
};

// Writes a trace of the current include file nesting to the 'out' stream.
void benchmark_include_file_trace( std::ostream& out, std::string fill) {
    std::size_t n = include_stack.size();
    for ( Include_stack::const_iterator i = include_stack.begin(); 
          i != include_stack.end(); ++i) 
    {
	out << std::setw(3) << n-- << ": '" << i->filename << "' line " 
            << i->linenumber << fill << std::endl;
    }
}

// This summy function is only here to suppress a warning for 
// an unused yyunput function 
void benchmark_parser_dummy_no_warn_() {
    yyunput( 1, 0);
}

// EOF

