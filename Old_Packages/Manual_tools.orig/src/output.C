/**************************************************************************
 
  output.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Output file management.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <output.h>
#include <fstream.h>
#include <list.h>
#include <html_config.h>
#include <html_error.h>
#include <macro_dictionary.h>
#include <string_conversion.h>
#include <cpp_formatting.h>

/* File and filename handling */
/* ========================== */
ostream* current_ostream  = 0;
string   current_filename = "<cout>";

ostream* pre_stream      = 0;
ostream* main_stream     = 0;
ostream* class_stream    = 0;
ostream* anchor_stream   = 0;
ostream* contents_stream = 0;
ostream* index_stream    = 0;

string   pre_main_filename;
string   main_filename = "<cout>";
string   class_filename;



/* Auxiliary functions for stream handling */
/* ======================================= */

bool exist_file( const string& name) {
    istream* in = new ifstream( name.c_str());
    if ( ! *in)
        return false;
    return true;
}

void assert_file_write( ostream& out, const string& name){
    if ( ! out) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot write to file `" 
	     << name << "'." << endl;
	exit(1);
    }
}

istream* open_file_for_read( const string& name){
    istream* in = new ifstream( name.c_str());
    if ( ! *in) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot open file `" << name
	     << "' for reading." << endl;
	exit(1);
    }
    return in;
}

istream* open_file_for_read_w_input_dirs( const string& name){

    if (name.at(0) == '/')  // an absolute path name is given
    {
       return open_file_for_read(name);
    }

    string::size_type first = 0;
    string::size_type last = 0;
    string dir = "";

    while (last < latex_conv_inputs.size())
    {
       last = latex_conv_inputs.find(':', first);
       if (last < latex_conv_inputs.size())
         dir = latex_conv_inputs.substr(first, last-first);
       else
         dir = latex_conv_inputs.substr(first, latex_conv_inputs.size()-first);
       assert_trailing_slash_in_path(dir);
       first = last+1;

       istream* in = new ifstream( (dir + name).c_str());
       if ( *in ) return in;
    }
    cerr << ' ' << endl 
	 << prog_name << ": error: cannot open file `" << name
	 << "' for reading." << endl;
    exit(1);
}

ostream* open_file_for_write( const string& name){
    ostream* out = new ofstream( name.c_str());
    if ( ! *out) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot open file `" << name
	     << "' for writing." << endl;
	exit(1);
    }
    return out;
}


// Output_file with stream operators
// -------------------------------------
Output_file current_output;

// Maintain a stack of output files
// --------------------------------
typedef list<Output_file> Output_file_stack;
Output_file_stack output_file_stack;

void push_current_output() {
    current_output = Output_file( current_ostream, current_filename);
    output_file_stack.push_front( current_output);
}

void pop_current_output() {
    if ( output_file_stack.empty())
	printErrorMessage( OutputStackEmptyError);
    else {
	current_output   = output_file_stack.front();
	current_ostream  = current_output.stream_ptr();
	current_filename = current_output.name();
	output_file_stack.pop_front();
    }
}

void set_current_output( const string& key) {
    if ( key == "premain")
	current_output = Output_file( pre_stream, pre_main_filename);
    else if ( key == "main")
	current_output = Output_file( main_stream, main_filename);
    else if ( key == "class")
	current_output = Output_file( class_stream, class_filename);
    else if ( key == "index")
	current_output = Output_file( index_stream,
				      macroX( "\\lciIndexFilename"));
    else if ( key == "toc")
	current_output = Output_file( contents_stream, 
				      macroX( "\\lciContentsFilename"));
    else if ( key == "anchor")
	current_output = Output_file( anchor_stream, 
				      macroX( "\\lciAnchorFilename"));
    else
	printErrorMessage( OutputStackKeyError);
    current_ostream  = current_output.stream_ptr();
    current_filename = current_output.name();
}

void push_current_output( const string& key) {
    push_current_output();
    set_current_output( key);
}


/* Filter a config file                    */
/* ======================================= */
/* Filtering means substituting of the variable names and */
/* skipping of braces if the variable is undefined.       */

void filter_config_file( istream& in, ostream& out) {
    char c;
    while( in.get(c)) {
        if ( c == '}')
	    continue;
        if ( c == '%' && in.get(c)) {
	    char d;
	    in.get(d);
	    if ( d == '{') {
	        bool is_empty = false;
		switch( c) {
		case 'c':
		    is_empty = class_name.empty();
		    break;
		case 'u':
		    is_empty = main_stream == &cout;
		    break;
		case 'd':
		    is_empty = macroX("\\lciManualDate").empty();
		    break;
		case 'r':
		    is_empty = macroX("\\lciManualRelease").empty();
		    break;
		case 't':
		    is_empty = macroX("\\lciManualTitle").empty();
		    break;
		case 'a':
		    is_empty = macroX("\\lciManualAuthor").empty();
		    break;
		default:
		    cerr << prog_name << ": warning: unknown variable scope %"
			 << c << "{...} in a config file found." << endl;
		    out.put('%');
		    out.put(c);
		    out.put(d);
		}
		if ( is_empty) {
		    // remove proper nested parantheses. Check for the
		    // escape sequence %{ or %}.
		    int nesting = 1; // count paranthesis nesting.
		    c = d;
		    while( in.get(d) && nesting > 0) {
		        if ( c != '%') {
			    if ( d == '{')
			        nesting ++;
			    if ( d == '}')
			        nesting --;
			}
		        c = d;
		    }
		}
	    } else {
		switch( c) {
		case '0':
		    out << prog_name;
		    break;
		case 'p':
		    out << prog_release;
		    break;
		case 'f':
		    out << current_filename;
		    break;
		case 'c':
		    if ( ! template_class_name.empty())
		        out << remove_font_commands( template_class_name);
		    break;
		case 'u':
		    if ( main_stream != &cout)
			if ( main_stream != pre_stream)
			    out <<  main_filename;
			else
			    out <<  pre_main_filename;
		    break;
		case 'd':
		    if ( ! macroX("\\lciManualDate").empty() ) 
		        out << macroX("\\lciManualDate");
		    break;
		case 'r':
		    if ( ! macroX("\\lciManualRelease").empty() ) 
		        out << macroX("\\lciManualRelease");
		    break;
		case 't':
		    if ( ! macroX("\\lciManualTitle").empty() ) 
		        out << macroX("\\lciManualTitle");
		    break;
		case 'a':
		    if ( ! macroX("\\lciManualAuthor").empty() ) 
		        out << macroX("\\lciManualAuthor");
		    break;
		case '%':
		case '{':
		case '}':
		    out.put(c);
		    break;
		default:
		    cerr << prog_name << ": warning: unknown variable %" << c 
			 << " in a config file found." << endl;
		    out.put('%');
		    out.put(c);
		}
		out.put(d);
	    }
	} else {
	    if ( in)
	        out.put(c);
	    else
		out.put('%');
	}
    }
}


/* Auxiliary functions to handle config files */
/* ========================================== */

istream* open_config_file( const string& name){
/*
    if ( config_switch == NO_SWITCH) {
        // check if the file exists in the current directory
        if ( exist_file( name))
	    return (open_file_for_read( name));
    }
*/
    return( open_file_for_read_w_input_dirs( config_path + name));
//    return( open_file_for_read_w_input_dirs(name));
}

void copy_and_filter_config_file( const string& name, ostream& out){
    istream* in  = open_config_file( name);
    filter_config_file( *in, out);
    delete in;
}

void copy_config_file( const string& name){
    istream* in  = open_config_file( name);
    ostream* out = open_file_for_write( tmp_path + name);
    filter_config_file( *in, *out);
    delete in;
    assert_file_write( *out, name);
    delete out;
}

void open_html( ostream& out) {
    istream* in = open_config_file( macroX( "\\lciManualHeader"));
    filter_config_file( *in, out);
    delete in;
}

void close_html( ostream& out) {
    istream* in = open_config_file( macroX( "\\lciManualFooter"));
    filter_config_file( *in, out);
    delete in;
}

// EOF //

