/**************************************************************************
 
  lex_include_impl.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Stack of files/strings to get scanned by a flex scanner.
              Internal implementation, only to be included by flex-file.
              Makes use of macros only defined in flex-file.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef LEX_INCLUDE_IMPL_H
#define LEX_INCLUDE_IMPL_H 1

#include <iomanip.h>
#include <macro_dictionary.h>

// Global variables used in scanner.
Include_stack       include_stack;
Include_stack_item* in_file       = 0;
Include_stack_item* in_string     = 0;


// Push current state. Init with new file and new_line_number.
void 
Include_stack::push_file( FILE* in, 
			  const string& name, 
			  size_t new_line_number) {
    if ( ! quiet_switch) {
	const string& cpath (macroX( "\\lciConfigPath"));
	if ( cpath.size() < name.size()) {
	    if ( name.substr(0, cpath.size()) == cpath)
		cerr << '[' << name.substr( cpath.size()) << flush;
	    else
		cerr << '[' << name << flush;
	} else
	    cerr << '[' << name << flush;
    }
    if( ! empty()) {
	m_stack.front().m_state = YY_START;
	m_stack.front().m_pbuf  = YY_CURRENT_BUFFER;
    }
    m_stack.push_front( Include_stack_item( in, name, new_line_number));
    in_file = in_string = &(m_stack.front());
    yyin = in;
    yy_switch_to_buffer( yy_create_buffer( yyin, YY_BUF_SIZE));
    insertInternalGlobalMacro( "\\lciInputFilename",    name);
    insertInternalGlobalMacro( "\\lciInputFilenameBase",basename_string(name));
    insertInternalGlobalMacro( "\\lciInputPath",        path_string( name));
}

// Push current state. Open and init with new file and new_line_number.
void 
Include_stack::push_file( const string& name, 
			  size_t new_line_number) {
    FILE* fin;
    if ( (fin = fopen( name.c_str(), "r")) == NULL) {
	cerr << ' ' << endl
	     << "Error: in `" << in_file->name() << "' line " 
	     << in_file->line() << ": cannot open include file `" 
	     << name << "' for reading." << endl;
    } else {
	push_file( fin, name, new_line_number);
    }
}

// Push current state. Open and init with new file plus optional
// extension string (suffix) and new_line_number.
void 
Include_stack::push_tex_file( const string& name, 
			      size_t new_line_number) {
    FILE* fin;
    if ( (fin = fopen( name.c_str(), "r")) != NULL)
	push_file( fin, name, new_line_number);
    else if ( (fin = fopen( (name + ".tex").c_str(), "r")) != NULL)
	push_file( fin, name, new_line_number);
    else if ( (fin = fopen( (name + ".sty").c_str(), "r")) != NULL)
	push_file( fin, name, new_line_number);
    else 
	cerr << ' ' << endl
	     << "Error: in `" << in_file->name() << "' line " 
	     << in_file->line() << ": cannot open file `" 
	     << name << "' for reading." << endl;
}

// Push current state. Init with new string and new_line_number.
void 
Include_stack::push_string( const string& name,
			    const string& s,
			    size_t new_line_number) {
    if ( s.empty()) {
	return;
    }
    if( ! empty()) {
	m_stack.front().m_state = YY_START;
	m_stack.front().m_pbuf  = YY_CURRENT_BUFFER;
    }
    m_stack.push_front( Include_stack_item( name, s, new_line_number));
    in_string = &(m_stack.front());
    yy_switch_to_buffer( yy_scan_string( (s + SEPARATOR).c_str()));
}

void 
Include_stack::pop() {
    assert( ! empty());
    if ( type() == Include_file) {
	if ( ! quiet_switch)
	    cerr << ']' << flush;
	fclose( file());
    }
    if ( size() > 1) {
	Stack::iterator i = m_stack.begin();
	++i;
	in_string = &*i;
	if ( type() == Include_file) {
	    in_file = in_string;
	    while ( in_file->type() != Include_file) {
		++i;
		if ( i == m_stack.end()) {
		    in_file = 0;
		    break;
		}
		in_file = &*i;
	    }
	    insertInternalGlobalMacro( "\\lciInputFilename",  in_file->name());
	    insertInternalGlobalMacro( "\\lciInputFilenameBase", 
				       basename_string( in_file->name()));
	    insertInternalGlobalMacro( "\\lciInputPath", 
				       path_string( in_file->name()));
	}
    } else {
	in_file = in_string = 0;
    }
    YY_BUFFER_STATE buf = YY_CURRENT_BUFFER;
    m_stack.pop_front();
    if ( ! empty())
	yy_switch_to_buffer( (YY_BUFFER_STATE)( m_stack.front().m_pbuf));
    yy_delete_buffer( buf);
    // BEGIN( state()); // Only in special cases by the caller.
}

void
Include_stack::trace( ostream& out) const {
    size_t n = m_stack.size();
    out << "Trace of the include-file stack:";
    if ( n == 0)
	out << " <empty>";
    out << endl;
    for ( Stack::const_iterator i = m_stack.begin(); i != m_stack.end(); ++i) {
	out << "    " << setw(3) << n--;
	if ( (*i).type() == Include_file)
	    out << "    File:  `";
	else
	    out << "  String:  `";
	out << (*i).name() << "' line " << (*i).line() << endl;
    }
}

#endif // LEX_INCLUDE_IMPL_H //
// EOF //
