/**************************************************************************
 
  lex_include.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Stack of files/strings to get scanned by a flex scanner.
              External declarations to be included in other files.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef LEX_INCLUDE_H
#define LEX_INCLUDE_H 1

#include <stdlib.h>
#include <stdio.h>
#include <mstring.h>
#include <list.h>

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
    typedef list<Include_stack_item>   Stack;
    typedef Stack::size_type  size_type;
    Stack   m_stack;

    // Push current state. Init with new file and new_line_number.
    bool                 push_file( FILE* in, 
				    const string& name,
				    size_t new_line_number = 1);

public:
    inline bool          empty() const { return m_stack.empty(); }
    inline size_type     size()  const { return m_stack.size(); }
    inline Include_type  type()  const { return m_stack.front().type(); }
    inline const string& name()  const { return m_stack.front().name(); }
    inline FILE*         file()  const { return m_stack.front().file(); }
    inline int           state() const { return m_stack.front().state(); }

    inline size_t        line()  const { return m_stack.front().line(); }
    inline size_t&       line()        { return m_stack.front().line(); }

    // Push current state. Open and init with new file and new_line_number.
    bool                 push_file( const string& name,
				    size_t new_line_number = 1);

    // Push current state. Open and init with new file plus optional
    // extension string (suffix) and new_line_number.
    bool                 push_tex_file( const string& name, 
					size_t new_line_number = 1);
    bool                 push_tex_file_w_input_dirs( const string& name, 
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
extern string              latex_conv_inputs; // defined in cc_extract_html.C

#endif // LEX_INCLUDE_H //
// EOF //
