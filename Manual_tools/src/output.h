/**************************************************************************

  output.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Output file management.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/

#ifndef OUTPUT_H
#define OUTPUT_H 1

#include <iostream>
#include <sstream>
#include <mstring.h>
#include <basic.h>

// Used in anchor rules to make anchor filter work for directory hierarchies
#define REPLACE_WITH_CURRENT_PATH_TOKEN "REPLACE_WITH_CURRENT_PATH"

using namespace std;

// Directory for the temporary files. A default is given.
// This directory is the output directory for this program. Usually,
// the script cc_manual_to_html makes a own tmp directory and removes
// it afterwards. The path must terminate with a slash.
extern string  tmp_path;


// Output_file with stream operators
// -------------------------------------

#define OUT(x) if ( m_out ) *m_out << x; return *this

class Output_file {
    ostream* m_out;
    string   m_name;
public:
    Output_file() : m_out(0), m_name( "<undef>") {}
    Output_file( ostream& out, const string& name)
	: m_out(&out), m_name(name) {}
    Output_file( ostream* out, const string& name)
	: m_out(out), m_name(name) {}
    bool           is_valid() const   { return m_out != 0;}
    ostream&       stream() const     { CC_Assert( m_out != 0); return *m_out;}
    ostream*       stream_ptr() const { return m_out;}
    const string&  name() const       { return m_name;}
    Output_file&  operator<<( char c)                   { OUT(c);}
    Output_file&  operator<<( const char* s)            { OUT(s);}
    Output_file&  operator<<( int a)                    { OUT(a);}
    Output_file&  operator<<( long l)                   { OUT(l);}
    Output_file&  operator<<( double d)                 { OUT(d);}
    Output_file&  operator<<( float f)                  { OUT(f);}
    Output_file&  operator<<( unsigned int a)           { OUT(a);}
    Output_file&  operator<<( unsigned long l)          { OUT(l);}
#ifdef _LONGLONG
    Output_file&  operator<<( long long l)              { OUT(l);}
    Output_file&  operator<<( unsigned long long l)     { OUT(l);}
#endif /* _LONGLONG */
    Output_file&  operator<<( void* p)                  { OUT(p);}
    Output_file&  operator<<( short i)                  { OUT(i);}
    Output_file&  operator<<( unsigned short i)         { OUT(i);}
    Output_file&  operator<<( ostream& (*f)(ostream&))  { OUT(f);}
    Output_file&  operator<<( ios& (*f)(ios&) )         { OUT(f);}

    operator const void*() const    { return m_out ? (void*)(*m_out) : 0; }
    Output_file&  flush() {
	if (m_out)
	    m_out->flush();
	return *this;
    }
    Output_file&  put(char c) {
	if (m_out)
	    m_out->put(c);
	return *this;
    }
    Output_file&  write(const char*  s,int n) {
	if (m_out)
	    m_out->write( s, n);
	return *this;
    }
};
#undef OUT

// actually, this is the output that _was_ the current output
// before calling push_current_output
extern Output_file current_output;

void push_current_output();
void pop_current_output();
void set_current_output( const string& key);
void push_current_output( const string& key);
void push_current_output_w_filename( const string& filename);

extern ostringstream* savebox_stream;

extern ostream* current_ostream;
extern string   current_filename;
extern string   current_basename;
extern string   current_rootname;
extern string   current_filepath;
extern string   current_uppath;

extern ostream* main_anchor_stream; // used in Chapter's
extern ostream* global_anchor_stream; // used for global files like TOC
extern ostream* anchor_stream; // the current, one of the above or class env.

extern ostream* pre_stream;
extern ostream* main_stream;
extern ostream* class_stream;
extern ostream* package_overview_stream;
extern ostream* contents_stream;
extern ostream* short_contents_stream;
extern ostream* comments_stream;
extern ostream* index_stream;
extern ostream* HREF_stream;
extern ostream* HREF_counter_stream;

extern string  pre_main_filename;
extern string  main_filename;
extern string  class_filename;

extern string  pre_main_basename;
extern string  main_basename;
extern string  class_basename;

extern string  pre_main_rootname;
extern string  main_rootname;
extern string  class_rootname;

extern string  pre_main_filepath;
extern string  main_filepath;
extern string  class_filepath;

extern string  pre_main_uppath;
extern string  main_uppath;
extern string  class_uppath;


/* Auxiliary functions for stream handling */
/* ======================================= */

bool     exist_file( const string& name); // also declared in input.h
void     assert_file_write( ostream& out, const string& name);
ostream* open_file_for_write( const string& name);
ostream* open_file_for_write_with_path( const string& name);
ostream* open_file_for_append( const string& name);
ostream* open_file_for_append_with_path( const string& name);
void     make_path( string path);

void     savestream_open   ( const string& name );
// returns NULL if not open
ostream* savestream_get    ( const string& name );
string   savestream_use    ( const string& name, const string& targetvarname );
void     savestream_close  ( const string& name );

void     minitoc_open();
void     minitoc_close();

#endif // OUTPUT_H //
// EOF //

