/**************************************************************************
 
  input.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Input file management and auxiliary functions.
              Migrated from 'lex_include.h':
              Stack of files/strings to get scanned by a flex scanner.
              External declarations to be included in other files.
              Implementation in 'lex_include_impl.h', to be included 
              only in the flex source file 'html_lex.yy'.
  System    : bison, flex, C++ (g++)
  Author    : (c) 2004 Lutz Kettner
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef INPUT_H
#define INPUT_H 1

#include <stdlib.h>
#include <stdio.h>
#include <mstring.h>
#include <list>
#include <iostream>
#include <basic.h>

using namespace std;


// Path for the HTML conversion tools for the default configuration files.
// This path will be compiled into the cc_extract_html program. It is set 
// in the Makefile. The same variable has to be configured in the 
// cc_manual_to_html script.
extern string  config_path;

// List of ':' separated paths (including '.' if needed) to check
// for LaTeX input files. Default is '.'. Can be set in the Makefile
// with the macro LATEX_CONV_INPUTS. Is overridden by the environment
// variable LATEX_CONV_INPUTS is it exists.
extern string  latex_conv_inputs;


/* Auxiliary functions for stream handling */
/* ======================================= */

// Returns true if file 'name' can be opened for read.
bool     exist_file( const string& name); // defined in output.cpp

// Checks if 'name' exists as readable file. Tries '.tex' and '.sty'
// suffixes in addition to the plain name. Returns found filename
// following the precedence order that TeX uses. Returns 'name' unchanged
// if file does not exist (has to be tested again separately).
string   find_filename_with_suffix( const string& name);

// Checks if 'name' exists as readable file, see find_filename_with_suffix.
// In addition, searches for 'name' exclusively in the subdirectories given
// in latex_conv_inputs separated by ':' if 'name' starts with a relative
// path. Returns "" if there is no file for that name.
string   find_filename_with_suffix_w_input_dirs( const string& name);

istream* open_file_for_read( const string& name);
istream* open_file_for_read_w_input_dirs( const string& name);
int      open_counter_file_for_read( const string& name);

// if a file 'name' should be included. in case of \includeonly{f}, 
// only 'f' is included. otherwise, all files are included
bool     is_to_be_included( const string& name );

// true iff there was some \includeonly statement with non-empty filename
bool     is_include_only();

// comma-separated list of files (excluding the .tex-filename extension)
void     include_only( const string& name );

/* Parts taken from lex_include.h          */
/* ======================================= */

enum Include_type { Include_file, Include_string};

class Include_stack;

class Include_stack_item {
    friend class Include_stack;

    Include_type    m_type;    // Include_file or Include_string;
    string          m_name;    // filename or name for string (macro etc.)
    int             m_state;   // Scanner state
    void*           m_pbuf;    // Opaque pointer to flex-buffer struct.
    // Used if type == Include_file
    FILE*           m_file;      
    size_t          m_line;
    // Used if type == Include_string
    string          m_str;
public:
    inline Include_type  type()  const { return m_type; }
    inline const string& name()  const { return m_name; }
    inline FILE*         file()  const { return m_file; }
    inline int           state() const { return m_state; }

    inline size_t        line()  const { return m_line; }
    inline size_t&       line()        { return m_line; }

    Include_stack_item() {}
    Include_stack_item( FILE* file, const string& name, size_t line)
	: m_type(Include_file),
	  m_name(name), 
	  m_state(0), 
	  m_pbuf(0), 
	  m_file(file), 
	  m_line(line)
    {}
    Include_stack_item( const string& name, const string& str, size_t line)
	: m_type(Include_string), 
	  m_name(name), 
	  m_state(0), 
	  m_pbuf(0), 
	  m_line(line), 
	  m_str(str)
    {}
};

class Include_stack {
    typedef std::list<Include_stack_item>   Stack;
    typedef Stack::size_type  size_type;
    Stack   m_stack;


public:
    inline bool          empty() const { return m_stack.empty(); }
    inline size_type     size()  const { return m_stack.size(); }
    inline Include_type  type()  const { return m_stack.front().type(); }
    inline const string& name()  const { return m_stack.front().name(); }
    inline FILE*         file()  const { return m_stack.front().file(); }
    inline int           state() const { return m_stack.front().state(); }

    inline size_t        line()  const { return m_stack.front().line(); }
    inline size_t&       line()        { return m_stack.front().line(); }

    // Push current state. Init with new file and new_line_number.
    bool                 push_file( FILE* in, 
				    const string& name,
				    size_t new_line_number = 1);
				    
    // Push current state. Open and init with new file and new_line_number.
    // Test for TeX file extensions and respect latex_conv_inputs.
    bool                 push_file( const string& name,
				    size_t new_line_number = 1);

    // Push current state. Init with new string and new_line_number.
    bool                 push_string( const string& name,
				      const string& s,
				      size_t new_line_number);

    void                 pop();
    void                 trace( ostream& out) const;

};

inline ostream&
operator<< ( ostream& out, const Include_stack& stack) {
    stack.trace( out);
    return out;
}

extern Include_stack       include_stack;
extern Include_stack_item* in_file;
extern Include_stack_item* in_string;


#endif // INPUT_H //
// EOF //

