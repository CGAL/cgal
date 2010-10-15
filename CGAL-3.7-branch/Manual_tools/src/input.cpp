/**************************************************************************
 
  input.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Input file management and auxiliary functions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 2004 Lutz Kettner
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <input.h>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <list>
#include <config.h>
#include <error.h>
#include <macro_dictionary.h>
#include <string_conversion.h>
#include <cpp_formatting.h>

using namespace std;



// Path for the HTML conversion tools for the default configuration files.
// This path will be compiled into the cc_extract_html program. It is set 
// in the Makefile. This string can be queried.
#ifndef LATEX_CONVERTER_CONFIG
#define LATEX_CONVERTER_CONFIG   ""
#endif
string config_path    = LATEX_CONVERTER_CONFIG;

// List of ':' separated paths (including '.' if needed) to check
// for LaTeX input files. Default is '.'. Can be set in the Makefile
// with the macro LATEX_CONVERTER_INPUTS. Is overridden by the environment
// variable LATEX_CONVERTER_INPUTS is it exists (at runtime)
#ifndef LATEX_CONVERTER_INPUTS
#define LATEX_CONVERTER_INPUTS   "."
#endif
string latex_conv_inputs   = LATEX_CONVERTER_INPUTS;

/* Auxiliary functions for stream handling */
/* ======================================= */

// Checks if 'name' exists as readable file. Tries '.tex' and '.sty'
// suffixes in addition to the plain name. Returns found filename
// following the precedence order that TeX uses. Returns 'name' unchanged
// if file does not exist (has to be tested again separately).
string find_filename_with_suffix( const string& name){
    FILE* fin = NULL;
    string found_name = name;
    string suffix = suffix_string( name);
    if (suffix == "tex" || suffix == "sty") {
        fin = fopen( name.c_str(), "r");
    } else if ( (fin = fopen( (name + ".tex").c_str(), "r")) != NULL) {
        found_name = name + ".tex";
    } else if ( (fin = fopen( (name + ".sty").c_str(), "r")) != NULL) {
        found_name = name + ".sty";
    }
    if ( fin != NULL)
        fclose(fin);
    return found_name;
}

// Checks if 'name' exists as readable file, see find_filename_with_suffix.
// In addition, searches for 'name' exclusively in the subdirectories given
// in latex_conv_inputs separated by ':' if 'name' starts with a relative
// path. Returns "" if there is no file for that name.
string find_filename_with_suffix_w_input_dirs( const string& name) {
    if ( name[0] == '/') { // an absolute path name is given
        string found_name = find_filename_with_suffix(name);
        if ( exist_file( found_name))
            return found_name;
        return "";
    }
    string::size_type first = 0;
    string::size_type last = 0;
    string dir = "";
    while (last < latex_conv_inputs.size()) {
        last = latex_conv_inputs.find(':', first);
        if (last < latex_conv_inputs.size())
            dir = latex_conv_inputs.substr(first, last-first);
        else
            dir = latex_conv_inputs.substr(first,
                                           latex_conv_inputs.size()-first);
        assert_trailing_slash_in_path(dir);
        first = last+1;
        string found_name = find_filename_with_suffix( dir + name);
        if ( exist_file( found_name))
            return found_name;
    }
    return "";
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
    string filename = find_filename_with_suffix_w_input_dirs( name);
    if ( filename == string("")) {
        cerr << ' ' << endl << "ERROR: "
             << prog_name << ": cannot open file `" << name
             << "' for reading." << endl;
        exit(1);
    }
    return open_file_for_read(filename);
}

int open_counter_file_for_read( const string& name){
    if (!exist_file(name.c_str())) {
       return 0;
    }
    istream* in = new ifstream( name.c_str());
    if ( ! *in) {
       return 0;
    }
    int i;
    *in >> i;
    delete in;
    return i;
}

static hash_set< string > files_to_be_included;


bool 
is_include_only() {
    return !files_to_be_included.empty();
}

bool     
is_to_be_included( const string& name ) {
    //std::cerr << "include? " << name << std::endl;
    return !is_include_only() || 
            files_to_be_included.find( name ) != files_to_be_included.end();
}


void     
include_only( const string& filename_list) {
    files_to_be_included.clear();
    
    string::size_type first = 0;
    string::size_type last = 0;
    string name = "";
    while( last < filename_list.size() ) {
        last = filename_list.find(',', first);
        if( last < filename_list.size() )
            name = filename_list.substr( first, last-first );
        else
            name = filename_list.substr( first,
                                        filename_list.size()-first );
        
        first = last+1;
        files_to_be_included.insert( name );
    }    
}



// EOF
