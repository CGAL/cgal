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
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef LEX_INCLUDE_IMPL_H
#define LEX_INCLUDE_IMPL_H 1

#include <iomanip.h>
#include <macro_dictionary.h>
#include <string_conversion.h>

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
    return true;
}

// Push current state. Open and init with new file and new_line_number.
bool
Include_stack::push_file( const string& name, 
			  size_t new_line_number) {
    FILE* fin;
    if ( (fin = fopen( name.c_str(), "r")) == NULL) {
	cerr << ' ' << endl 
	     << "*** Error: cannot open file `" << name << "' for reading.";
	printErrorMessage( FileReadOpenError);
	return false;
    }
    return push_file( fin, name, new_line_number);
}

// Push current state. Open and init with new file plus optional
// extension string (suffix) and new_line_number.
bool
Include_stack::push_tex_file( const string& name, 
			      size_t new_line_number) {
    FILE* fin;
    string suffix = suffix_string( name);
    if ( (suffix == "tex" || suffix == "sty") &&
	 ( (fin = fopen( name.c_str(), "r")) != NULL))
	return push_file( fin, name, new_line_number);
    if ( (fin = fopen( (name + ".tex").c_str(), "r")) != NULL)
	return push_file( fin, name + ".tex", new_line_number);
    if ( (fin = fopen( (name + ".sty").c_str(), "r")) != NULL)
	return push_file( fin, name + ".sty", new_line_number);
    if ( (fin = fopen( name.c_str(), "r")) != NULL)
	return push_file( fin, name, new_line_number);
    return false;
}



bool
Include_stack::push_tex_file_w_input_dirs( const string& name, 
			                  size_t new_line_number)
{
    string::size_type first = 0;
    string::size_type last = 0;
    string dir;

    // check for absolute path 
    if (name[0] == '/') 
      if (push_tex_file(name, new_line_number))
      {
        return true;
      }
      else
      {
         cerr << ' ' << endl 
	      << "*** Error: cannot open file `" << name << "' for reading.";
         printErrorMessage( FileReadOpenError);
         return false;
      }


    while (last < latex_conv_inputs.size())
    {
       last = latex_conv_inputs.find(':', first);
       if (last < latex_conv_inputs.size())
          dir = latex_conv_inputs.substr(first, last-first);
       else
          dir = latex_conv_inputs.substr(first, latex_conv_inputs.size()-first);
       assert_trailing_slash_in_path(dir);
       first = last+1;
       if (push_tex_file(dir+name, new_line_number)) return true;
    }
    cerr << ' ' << endl 
	 << "*** Error: cannot open file `" << name << "' for reading.";
    printErrorMessage( FileReadOpenError);
    return false;
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
