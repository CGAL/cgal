/**************************************************************************
 
  lex_include_impl.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Stack of files/strings to get scanned by a flex scanner.
              Internal implementation, only to be included by flex-file.
              Makes use of macros only defined in flex-file.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef LEX_INCLUDE_IMPL_H
#define LEX_INCLUDE_IMPL_H 1

#include <iomanip>
#include <macro_dictionary.h>
#include <string_conversion.h>

using namespace std;

// Global variables used in scanner.
Include_stack       include_stack;
Include_stack_item* in_file       = 0;
Include_stack_item* in_string     = 0;


// Push current state. Init with new file and new_line_number.
bool
Include_stack::push_file( FILE* in, 
			  const string& name, 
			  size_t new_line_number) {
    if ( ! quiet_switch) {
	const string& cpath (macroX( "\\lciConfigPath"));
	if ( cpath.size() < name.size() && 
             name.substr(0, cpath.size()) == cpath )
            cerr << '[' << name.substr( cpath.size()) << flush;
	else
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
    insertInternalGlobalMacro( "\\lciInputFilename", name);
    insertInternalGlobalMacro( "\\lciInputBasename", basename_string(name));
    insertInternalGlobalMacro( "\\lciInputRootname", rootname_string(name));
    string path = path_string( name);
    insertInternalGlobalMacro( "\\lciInputPath",     path);
    insertInternalGlobalMacro( "\\lciInputUppath",   uppath_string( path));
    return true;
}

// Push current state. Open and init with new file and new_line_number.
bool
Include_stack::push_file( const string& name, 
			  size_t new_line_number) {
    string filename = find_filename_with_suffix_w_input_dirs( name);
    if ( filename == string("")) {
        cerr << ' ' << endl << "ERROR: "
             << prog_name << ": cannot open file `" << name
             << "' for reading." << endl;
	printErrorMessage( FileReadOpenError);
	return false;
    }
    FILE* fin = fopen( filename.c_str(), "r");
    CC_Assert( fin != NULL);
    return push_file( fin, filename, new_line_number);
}


// Push current state. Init with new string and new_line_number.
bool
Include_stack::push_string( const string& name,
			    const string& s,
			    size_t new_line_number) {
    if ( s.empty()) {
	return true;
    }
    if( ! empty()) {
	m_stack.front().m_state = YY_START;
	m_stack.front().m_pbuf  = YY_CURRENT_BUFFER;
    }
    m_stack.push_front( Include_stack_item( name, s, new_line_number));
    in_string = &(m_stack.front());
    yy_switch_to_buffer( yy_scan_string( (s + SEPARATOR).c_str()));
    return true;
}

void 
Include_stack::pop() {
    if ( empty()) {
	printErrorMessage( IncludeStackUnderflowError);
	return;
    }
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
	    insertInternalGlobalMacro( "\\lciInputBasename", 
				       basename_string( in_file->name()));
	    insertInternalGlobalMacro( "\\lciInputRootname", 
				       rootname_string( in_file->name()));
            string path = path_string( in_file->name());
            insertInternalGlobalMacro( "\\lciInputPath",   path);
            insertInternalGlobalMacro( "\\lciInputUppath",uppath_string(path));
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
