/**************************************************************************
 
  output.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Output file management.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef OUTPUT_H
#define OUTPUT_H 1

#include <iostream.h>
#include <mstring.h>
#include <basic.h>

#define OUT(x) if ( m_out) *m_out << x; return *this

// Output_file with stream operators
// -------------------------------------
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

extern Output_file current_output;

void push_current_output();
void pop_current_output();
void set_current_output( const string& key);
void push_current_output( const string& key);

extern ostream* current_ostream;
extern string   current_filename;

extern ostream* pre_stream;
extern ostream* main_stream;
extern ostream* class_stream;
extern ostream* anchor_stream;
extern ostream* contents_stream;
extern ostream* index_stream;

extern string  pre_main_filename;
extern string  main_filename;
extern string  class_filename;

extern string  config_path;        // defined in cc_extract_html.C
extern string  latex_conv_inputs;  // 
extern string  tmp_path;           //

/* Auxiliary functions for stream handling */
/* ======================================= */

bool     exist_file( const string& name);
void     assert_file_write( ostream& out, const string& name);
istream* open_file_for_read( const string& name);
istream* open_file_for_read_w_input_dirs( const string& name);
ostream* open_file_for_write( const string& name);


/* Filter a config file                    */
/* ======================================= */
/* Filtering means substituting of the variable names and */
/* skipping of braces if the variable is undefined.       */

void filter_config_file( istream& in, ostream& out);

/* Auxiliary functions to handle config files */
/* ========================================== */

istream* open_config_file( const string& name);
void     copy_and_filter_config_file( const string& name, ostream& out);
void     copy_config_file( const string& name);

void open_html( ostream& out);
void close_html( ostream& out);



#endif // OUTPUT_H //
// EOF //

