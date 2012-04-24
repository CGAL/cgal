/**************************************************************************

  output.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Output file management.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/

#include <output.h>
#include <sstream>
#include <fstream>
#include <cassert>
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

// Directory for the temporary files. A default is given.
// This directory is the output directory for this program. Usually,
// the script cc_manual_to_html makes a own tmp directory and removes
// it afterwards. The path must terminate with a slash.
string tmp_path       = "/usr/tmp/";

/* File and filename handling */
/* ========================== */
ostream* current_ostream  = 0;
string   current_filename = "<cout>";
string   current_basename = current_filename;
string   current_rootname = current_filename;
string   current_filepath;
string   current_uppath;

ostream* main_anchor_stream    = 0; // used in Chapter's
ostream* global_anchor_stream  = 0; // used for global files like TOC
ostream* anchor_stream         = 0; // the current, one of the above or class

ostream* pre_stream          = 0;
ostream* main_stream         = 0;
ostream* package_overview_stream  = 0;
ostream* class_stream        = 0;
ostream* contents_stream     = 0;
ostream* short_contents_stream = 0;
ostream* comments_stream = 0;
ostream* index_stream = 0;
ostream* HREF_stream = 0;
ostream* HREF_counter_stream = 0;
ostream* minitoc_stream = 0;

ostringstream* savebox_stream = new ostringstream;

string   pre_main_filename;
string   main_filename = "<cout>";
string   class_filename;

string   pre_main_basename;
string   main_basename = main_filename;
string   class_basename;

string   pre_main_rootname;
string   main_rootname = main_filename;
string   class_rootname;

string   pre_main_filepath;
string   main_filepath;
string   class_filepath;

string   pre_main_uppath;
string   main_uppath;
string   class_uppath;

/* Auxiliary functions for stream handling */
/* ======================================= */

bool exist_file( const string& name) {
    ifstream in( name.c_str());
    if ( ! in)
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

void make_path( string path){
    if ( path.size() < 2)
        return;
    path.replace( path.size() - 1, 1, "");
    struct stat st;
    int s = stat( path.c_str(), &st);
    if ( s != 0) {
        make_path( path_string( path));
        if ( mkdir( path.c_str(), 0755)) {
            cerr << ' ' << endl
                 << prog_name << ": error: cannot create directory `" << path
                 << "'." << endl;
            exit(1);
        }
    }
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

ostream* open_file_for_write_with_path( const string& name){
    make_path( path_string( name));
    return open_file_for_write( name);
}

ostream* open_file_for_append( const string& name){
    ostream* out = new ofstream( name.c_str(), std::ios::out | std::ios::app);
    if ( ! *out) {
        cerr << ' ' << endl
	     << prog_name << ": error: cannot open file `" << name
	     << "' for writing." << endl;
	exit(1);
    }
    return out;
}

ostream* open_file_for_append_with_path( const string& name){
    make_path( path_string( name));
    return open_file_for_append( name);
}


// Output_file with stream operators
// -------------------------------------
Output_file current_output;

// Maintain a stack of output and anchor files
// -------------------------------------------
typedef std::list<Output_file> Output_file_stack;
typedef std::list<ostream*>    Anchor_file_stack;

Output_file_stack output_file_stack;
Anchor_file_stack anchor_file_stack;

void push_current_output() {
    if (contents_stream)
        (*contents_stream) << flush;
    current_output = Output_file( current_ostream, current_filename);
    output_file_stack.push_front( current_output );
    anchor_file_stack.push_front( anchor_stream );
}

void pop_current_output() {
    encode_backslash( true );
    if( output_file_stack.empty() )
        printErrorMessage( OutputStackEmptyError);
    else {
        if( contents_stream )
            (*contents_stream) << flush;
        current_output   = output_file_stack.front();
        anchor_stream    = anchor_file_stack.front( );
        current_ostream  = current_output.stream_ptr();
        current_filename = current_output.name();
        current_basename = basename_string( current_filename );
        current_rootname = rootname_string( current_basename );
        current_filepath = path_string( current_filename );
        current_uppath   = uppath_string( current_filepath );
        insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
        insertInternalGlobalMacro( "\\lciOutputBasename",current_basename);
        insertInternalGlobalMacro( "\\lciOutputRootname",current_rootname);
        insertInternalGlobalMacro( "\\lciOutputPath",    current_filepath);
        insertInternalGlobalMacro( "\\lciOutputUppath",  current_uppath);
        output_file_stack.pop_front();
        anchor_file_stack.pop_front();
    }
}

void set_current_output( const string& key) {
    bool new_filename = true;
    if ( key == "premain") {
        current_output = Output_file( pre_stream, pre_main_filename);
        anchor_stream = global_anchor_stream;
    } else if ( key == "main") {
        current_output = Output_file( main_stream, main_filename);
        anchor_stream = main_anchor_stream;
    } else if ( key == "class") {
        current_output = Output_file( class_stream, class_filename);
    } else if ( key == "packages") {
        current_output = Output_file( package_overview_stream,
                                      macroX( "\\lciPkgOverviewFilename") );
        anchor_stream = global_anchor_stream;
    } else if ( key == "toc") {
        current_output = Output_file( contents_stream,
                                      macroX( "\\lciContentsFilename") );
        anchor_stream = global_anchor_stream;
    } else if ( key == "comments") {
        current_output = Output_file( comments_stream, string("comments.xml") );
        new_filename = false;
        encode_backslash( false );
    } else if ( key == "index") {
        current_output = Output_file( index_stream,
                                      macroX( "\\lciIndexFilename") );
        anchor_stream = global_anchor_stream;
        new_filename = false;
    } else if ( key == "anchor") {
        current_output = Output_file( anchor_stream,
                                      macroX( "\\lciAnchorFilename") );
        new_filename = false;
    } else if ( key == "savebox" ) {
        current_output = Output_file( savebox_stream, "savebox" );
        new_filename = false;
    } else if ( key == "savestream" ) {
        string savestream_name = macroX( "\\lciSaveStreamName" );
        ostream* savestream = savestream_get( savestream_name );
        if( savestream != NULL )
          current_output = Output_file( savestream, savestream_name );
        else {
          std::cerr << "!! Error: savestream \"" << savestream_name << "\" unknown!" << std::endl;
          return;
        }
        //new_filename = false;
        anchor_stream = global_anchor_stream;
    } else if ( key == "minitoc" ) {
        string filename = macroX( "\\lciMinitocFilename" );
        current_output = Output_file( minitoc_stream, filename );
        new_filename = false;
    }
    else {
        printErrorMessage( OutputStackKeyError);
    }
    current_ostream  = current_output.stream_ptr();
    if ( new_filename) {
        current_filename = current_output.name();
        current_basename = basename_string( current_filename);
        current_rootname = rootname_string( current_basename);
        current_filepath = path_string( current_filename);
        current_uppath   = uppath_string( current_filepath);
        insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
        insertInternalGlobalMacro( "\\lciOutputBasename",current_basename);
        insertInternalGlobalMacro( "\\lciOutputRootname",current_rootname);
        insertInternalGlobalMacro( "\\lciOutputPath",    current_filepath);
        insertInternalGlobalMacro( "\\lciOutputUppath",  current_uppath);
    }
}

void push_current_output( const string& key) {
    push_current_output();
    set_current_output( key);
}

void push_current_output_w_filename( string filename,
                                     const string path) {
    push_current_output();
    current_filename = filename;
    if( !path.empty() ) {
      filename = path + filename;
    }
    current_ostream  = open_file_for_write( filename);
    current_basename = basename_string( filename);
    current_rootname = rootname_string( current_basename);
    current_filepath = path_string( filename);
    current_uppath   = uppath_string( current_filepath);
    insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
    insertInternalGlobalMacro( "\\lciOutputBasename",current_basename);
    insertInternalGlobalMacro( "\\lciOutputRootname",current_rootname);
    insertInternalGlobalMacro( "\\lciOutputPath",    current_filepath);
    insertInternalGlobalMacro( "\\lciOutputUppath",  current_uppath);
}


typedef hash_map< string, ostream* > Savestream_table;
Savestream_table savestream_table;

void
savestream_open( const string& name ) {
  Savestream_table::iterator it = savestream_table.find( name );
  if( it == savestream_table.end() || savestream_table[ name ] == NULL)
    savestream_table[ name ] = new ostringstream;
  else
    std::cerr << "!! Error: savestream \"" << name << "\" already open!" << std::endl;
}

ostream*
savestream_get( const string& name ) {
  Savestream_table::iterator it = savestream_table.find( name );
  if( it != savestream_table.end() )
    return it->second;
  else
    return NULL;
}

string
savestream_use( const string& name, const string& targetvarname ) {
  ostream *out = savestream_get( name );
  if( out != NULL ) {
    ostringstream *sout = dynamic_cast<ostringstream*>(out);
    assert( sout != NULL );
    return "\\gdef\\" + targetvarname + "{" + sout->str() + "}";
  } else {
    std::cerr << "!! Warning: savestream '" << name << "' not found" << std::endl;
  }
  return string();
}

void
savestream_close( const string& name ) {
  ostream *out = savestream_get( name );
  if( out != NULL ) {
    delete out;
    savestream_table[name] = NULL;
  }
}

void
minitoc_open() {
   if( minitoc_stream != NULL ) {
     std::cerr << "!! Error: minitoc opened twice ... current filepath: " << current_filepath << std::endl;
     delete minitoc_stream;
   }
   string filename = current_filename + ".minitoc";
   insertInternalGlobalMacro( "\\lciMinitocFilename",filename);
   
   minitoc_stream = open_file_for_write_with_path( tmp_path + filename );
}

void
minitoc_close() {
   if( minitoc_stream != NULL ) {
     delete minitoc_stream;
     minitoc_stream = NULL;
   }
}

// EOF //

