/**************************************************************************
 
  cc_extract_html.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : main program, command line parameter parsing
  System    : bison, flex, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stream.h>
#include <fstream.h>
#include <strstream.h>
#include <ctype.h>

#include <mstring.h>

#include <lex_include.h>
#include <macro_dictionary.h>
#include <internal_macros.h>
#include <buffer.h>
#include <html_config.h>
#include <html_syntax.h>
#include <html_error.h>
#include <string_conversion.h>
#include <output.h>


/* Program name and release          */
/* ================================= */
const string prog_name    = "cc_extract_html";
const string prog_release = "$Revision$";

// Path for the HTML conversion tools for the default configuration files.
// This path will be compiled into the cc_extract_html program. It is set 
// in the Makefile. The same variable has to be configured in the 
// cc_manual_to_html script.
#ifndef LATEX_CONVERTER_CONFIG
#define LATEX_CONVERTER_CONFIG   ""
#endif
#ifndef LATEX_CONV_INPUTS
#define LATEX_CONV_INPUTS   ""
#endif
string config_path    = LATEX_CONVERTER_CONFIG;
string latex_conv_inputs   = LATEX_CONV_INPUTS;


// Directory for the temporary files. A default is given.
// This directory is the output directory for this program. Usually, 
// the script cc_manual_to_html makes a own tmp directory and removes 
// it afterwards. The path must terminate with a slash.
string tmp_path       = "/usr/tmp/";

/* Constant string used for referencing. */
/* ===================================== */
/* This constant must be doubly quoted since it is subject of another */
/* C compiler pass during generation of the link generator.           */
const string reference_icon = "<IMG SRC=\\\"cc_ref_up_arrow.gif\\\" "
            "ALT=\\\"reference\\\" WIDTH=\\\"10\\\" HEIGHT=\\\"10\\\">";

/* Configurable command line options */
/* ================================= */
Switch  V_switch           = NO_SWITCH;

Switch  trace_switch       = NO_SWITCH;
Switch  line_switch        = NO_SWITCH;

Switch  config_switch      = NO_SWITCH;
Switch  quiet_switch       = NO_SWITCH;
Switch  macro_def_switch   = NO_SWITCH;
Switch  macro_exp_switch   = NO_SWITCH;
Switch  macro_def2_switch  = NO_SWITCH;
Switch  macro_exp2_switch  = NO_SWITCH;
Switch  sty_macro_switch   = NO_SWITCH;
Switch  stack_trace_switch = NO_SWITCH;

Switch  noheader_switch    = NO_SWITCH;
Switch  onlyheader_switch  = NO_SWITCH;


/* Config filename:                  */
/* ================================= */
string sty_filename      = "default.sty";


/* Configurable command line options */
/* ================================= */

// Manual date, release, title, and author used in filtering the config files.
string manual_date;
string manual_release;
string manual_title;
string manual_author;


void init_commandline_args() {
    insertInternalGlobalMacro( "\\lciConfigPath",     config_path);
    insertInternalGlobalMacro( "\\lciTmpPath",        tmp_path);
    insertInternalGlobalMacro( "\\lciReferenceIcon",  reference_icon);
    insertInternalGlobalMacro( "\\lciManualDate",     manual_date);
    insertInternalGlobalMacro( "\\lciManualRelease",  manual_release);
    insertInternalGlobalMacro( "\\lciManualTitle",    manual_title);
    insertInternalGlobalMacro( "\\lciManualAuthor",   manual_author);
    insertInternalGlobalMacro( "\\lciIfVersionFlag",  V_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfQuietFlag",    quiet_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfMacroDefFlag", macro_def2_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfMacroExpFlag", macro_exp2_switch 
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfNoHeaderFlag", noheader_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciOutputFilename", "<cout>");
    insertInternalGlobalMacro( "\\lciMainFilename",   "<cout>");
}




int text_block_length( const Buffer_list& T) {
    int l = 0;
    for ( Buffer_const_iterator words = T.begin(); words != T.end(); ++words) {
        if ( (*words)->is_space())
	    ++l;
	else
	    l += (*words)->size() - 1;
    }
    return l;
}

char* text_block_to_string( const Buffer_list& T) {
    char* string = new char[ text_block_length(T) + 1];
    string[0] = '\0';
    for ( Buffer_const_iterator words = T.begin(); words != T.end(); ++words) {
        if ( (*words)->is_space())
	    strcat( string, " ");
	else
	    strcat( string, (*words)->begin());
    }
    return string;
}

bool is_text_block_empty( const Buffer_list& T) {
    for ( Buffer_const_iterator words = T.begin(); words != T.end(); ++words) {
        if ( ! (*words)->is_space())
	    return false;
    }
    return true;
}

void print_html_text_block( ostream &out, const Buffer_list& T) {
    for ( Buffer_const_iterator words = T.begin(); words != T.end(); ++words) {
	const char* s = (*words)->begin();
	while ( *s) {
	    if ( *s != SEPARATOR)
		out << *s;
	    ++s;
	}
    }
}



/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleText( const Buffer_list& T) {
    if ( ! current_ostream)
	return;
    print_html_text_block( *current_ostream, T);
}

void handleBuffer( const Buffer& TT) {
    if ( current_ostream)
	*current_ostream << TT.begin();
}

void handleString( const char* s) {
    if ( current_ostream)
	*current_ostream << s;
}

void handleString( const string& s) {
    if ( current_ostream)
	*current_ostream << s;
}

void handleChar( char c) {
    if ( current_ostream)
	*current_ostream << c;
}


/* main */
/* ==== */
 
#define MaxParameters           1000
#define MaxOptionalParameters    999
#define ErrParameters          10000
 
/* this macro opens a block, in which the switch is detected */
/* it must be closed with the macro endDetect()              */
#define detectSwitch( var, text) \
    if ( (( argv[i][0] == '/' ) || ( argv[i][0] == '-' ) || \
          ( argv[i][0] == '+' )) && ( strcmp( text, argv[i]+1) == 0)) { \
        if ( argv[i][0] == '+' ) \
            var = PLUS_SWITCH; \
        else \
            var = MINUS_SWITCH;
 
#define endDetect() \
        if ( nParameters <= MaxParameters ) \
            continue; \
        else \
            break; \
    }
 
 
 
/* >main: main function with standard unix parameter input */
/* ------------------------------------------------------- */
 
main( int argc, char **argv) {
    // Check environment:
    char* s = getenv("LATEX_CONV_CONFIG");
    if ( s)
	config_path = s;
    s = getenv("LATEX_CONV_INPUTS");
    if ( s)
       latex_conv_inputs = s;
    else
       latex_conv_inputs = "."; // if the environment variable is not set,
                                // try to get all files from current directory

    int i;
    int nParameters = 0;
    typedef const char* ParamType;
    ParamType parameters[ MaxParameters + 1];

    Switch help_switch  = NO_SWITCH;
    Switch dummy_switch;
    for (i = 1; i < argc; i++) {

        /* check switches */
        detectSwitch( dummy_switch, "date");
            i++;
            if ( i < argc)
	        manual_date = argv[i];
	    else {
	        cerr << "*** Error: option -date needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "release");
            i++;
            if ( i < argc)
	        manual_release = argv[i];
	    else {
	        cerr << "*** Error: option -release needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "title");
            i++;
            if ( i < argc)
	        manual_title = argv[i];
	    else {
	        cerr << "*** Error: option -title needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "author");
            i++;
            if ( i < argc)
	        manual_author = argv[i];
	    else {
	        cerr << "*** Error: option -author needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( config_switch, "config");
            i++;
            if ( i < argc) {
	        config_path = argv[i];
	    } else {
	        cerr << "*** Error: option -config needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "sty");
            i++;
            if ( i < argc) {
	        sty_filename = argv[i];
	    } else {
	        cerr << "*** Error: option -sty needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "tmp");
            i++;
            if ( i < argc) {
	        tmp_path = argv[i];
	    } else {
	        cerr << "*** Error: option -tmp needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "header");
            i++;
            if ( i < argc) {
		string s = argv[i];
		assert_trailing_slash_in_path( s);
		insertGlobalMacro( "\\lciHeaderPath", 
				   "<command line option>", 0, s);
	    } else {
	        cerr << "*** Error: option -header needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
	// The following block remains for compatibility.
        detectSwitch( dummy_switch, "cgal_dir");
            i++;
            if ( i < argc) {
		string s = argv[i];
		assert_trailing_slash_in_path( s);
		insertGlobalMacro( "\\lciHeaderPath", 
				   "<command line option>", 0, s);
	    } else {
	        cerr << "*** Error: option -cgal_dir needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "main");
            i++;
            if ( i < argc) {
	        pre_main_filename = argv[i];
	    } else {
	        cerr << "*** Error: option -main needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( quiet_switch, "quiet");
        endDetect();
        detectSwitch( trace_switch, "trace");
        endDetect();
        detectSwitch( line_switch, "line");
        endDetect();
        detectSwitch( macro_def2_switch, "macrodef");
        endDetect();
        detectSwitch( macro_exp2_switch, "macroexp");
        endDetect();
        detectSwitch( sty_macro_switch, "stymacro");
        endDetect();
        detectSwitch( stack_trace_switch, "stacktrace");
        endDetect();
        detectSwitch( noheader_switch, "noheader");
        endDetect();
        detectSwitch( onlyheader_switch, "onlyheader");
        endDetect();

        detectSwitch( V_switch, "V");
	    cout << prog_name << " " << prog_release << " (c) Lutz Kettner" 
		 << endl;
	    cout << "Using: ";
        endDetect();
        detectSwitch( help_switch, "h");
        endDetect();
        detectSwitch( help_switch, "H");
        endDetect();
        detectSwitch( help_switch, "help");
        endDetect();

	// check for unknown command line option
	if ( argv[i][0] == '-' ) {
	    cerr << "*** Error: unknown command line option `" << argv[i] 
		 << "'." << endl;
	    nParameters = ErrParameters;
	}

        /* else get standard or optional parameters */
        if ( nParameters < MaxParameters ) {
            parameters[nParameters ++] = argv[i];
            continue;
        }
 
        nParameters = ErrParameters;
        break;
    }
    (void)(dummy_switch);  // simulate a use of 'dummy_switch'.

    if (! V_switch && 
        ((nParameters < MaxParameters - MaxOptionalParameters) ||
         (nParameters > MaxParameters) || (help_switch != NO_SWITCH))) {
        if (help_switch == NO_SWITCH)
            cerr << "*** Error: in parameter list" << endl;
	cerr << prog_name << " " << prog_release << " (c) Lutz Kettner" 
	     << endl;
        cerr << "Usage: " << prog_name << " [<options>] <TeX-files...>" 
	     << endl;
        cerr << "       -V                  prints version message." <<endl;
        cerr << "       -date     <text>    set a date for the manual." <<endl;
        cerr << "       -release  <text>    set a release number for the "
                                           "manual." << endl;
        cerr << "       -title    <text>    set a title text for the manual." 
	     << endl;
        cerr << "       -author   <text>    set an author address (email) for "
	                                   "the manual." << endl;
        cerr << "       -config   <dir>     set the path where to find the "
                                           "config files." << endl;
        cerr << "       -tmp      <dir>     set the path where to put the "
                                           "output files." << endl;
        cerr << "       -header   <dir>     set the path to the C "
                                           "header files." << endl;
        cerr << "       -main     <file>    main filename for the part before"
                                           " any chapter." << endl;
        cerr << "       -sty      <style>   use style file." << endl;
        cerr << "       -quiet              no output, no warnings for "
                                            "unknown macros." << endl;
        cerr << "       -macrodef           trace macro definitions." << endl;
        cerr << "       -macroexp           trace macro expansions." << endl;
        cerr << "       -stymacro           trace style macros as well."<<endl;
        cerr << "       -stacktrace         stack trace for each error."<<endl;
        cerr << "       -noheader           no HTML header for contents.html "
                "and index." << endl;
        cerr << "       -onlyheader         convert config files instead of "
                "TeX-files." << endl;
        cerr << "       -trace              set `yydebug' for bison to true."
             << endl;
        cerr << "       -line               echo currently parsed line "
                                           "numbers to cerr." << endl;
        cerr << "(TeX-files with suffix:  .tex  or  .bbl  possible.)"  
	     << endl;
        exit(1);
    }

    // Prepare proper format of path arguments.
    assert_trailing_slash_in_path( tmp_path);
    assert_trailing_slash_in_path( config_path);
    config_path += "html/";
 
    // Initialization
    if ( ! quiet_switch && ! V_switch)
	cerr << '[' << prog_name << ' ' << prog_release << endl;
    if ( sty_macro_switch) {
	macro_def_switch = macro_def2_switch;
	macro_exp_switch = macro_exp2_switch;
	if ( trace_switch)
	    yydebug = 1;
    }
    init_commandline_args();
    init_internal_macros();
    current_ostream  = 0;

    if (!include_stack.push_tex_file_w_input_dirs(config_path + sty_filename))
	exit(1);

    yyparse();
    if (V_switch) {
	cerr << endl;
	exit(0);
    }

    main_stream   = &cout;

    if ( onlyheader_switch) {
	for ( i = 0; i < nParameters; i++)
	    copy_config_file( parameters[i]);
	index_stream = open_file_for_write(tmp_path + 
					   macroX( "\\lciIndexFilename"));
	write_headers_to_index( *index_stream);
	assert_file_write( *index_stream, 
			   macroX( "\\lciIndexFilename"));
	delete index_stream;
	if ( ! quiet_switch)
	    cerr << ']' << endl;
	return 0;
    }

    if ( ! noheader_switch) {
	// Filter config files for the index.
	copy_config_file( macroX( "\\lciIndexHeader"));
	copy_config_file( macroX( "\\lciIndexFooter"));
    }

    // Prepare several streams:
    // anchor_stream = new ofstream( anchor_filename, ios::out | ios::app);
    anchor_stream   = open_file_for_write( tmp_path +
					   macroX( "\\lciAnchorFilename"));
    contents_stream = open_file_for_write( tmp_path +
					   macroX( "\\lciContentsFilename"));
/*
    if ( ! noheader_switch)
	copy_and_filter_config_file( macroX( "\\lciTocHeader"), 
				     *contents_stream);
*/

    index_stream    = open_file_for_write( tmp_path +
					   macroX( "\\lciIndexFilename"));
    if ( ! noheader_switch)
	write_headers_to_index( *index_stream);

    if ( ! pre_main_filename.empty()) {
	pre_stream = open_file_for_write( tmp_path + pre_main_filename);
	current_filename = pre_main_filename;
	open_html( *pre_stream);
    }
    
    macro_def_switch = macro_def2_switch;
    macro_exp_switch = macro_exp2_switch;
    if ( trace_switch)
	yydebug = 1;

    for ( i = 0; i < nParameters; i++) {
	if ( ! pre_main_filename.empty()) {
	    main_stream      = pre_stream;
	    main_filename    = pre_main_filename;
	    current_ostream  = main_stream;
	    current_filename = pre_main_filename;
	    insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
	    insertInternalGlobalMacro( "\\lciMainFilename",  main_filename);
	} else {
	    main_stream      = &cout;
	    current_ostream  = main_stream;
	    current_filename = main_filename;
	    insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
	}
	if ( include_stack.push_tex_file_w_input_dirs( parameters[i]))
	    yyparse();

	include_stack.push_string( "<end of conversion>", 
				   "\\lciCheckNestingScopes", 
				   0);
	include_stack.push_string( "<end of conversion>", 
				   "\\lciEndOfConversion", 
				   0);
	yyparse();

	assert_file_write( *main_stream, main_filename);
	if ( main_stream != &cout && main_stream != pre_stream) {
	    close_html( *main_stream);
	    assert_file_write( *main_stream, main_filename);
	    delete   main_stream;
	    main_stream = &cout;
	    main_filename = "<cout>";
	}
    }
    if ( ! quiet_switch)
	cerr << ']' << endl;

    if ( ! pre_main_filename.empty()) {
        close_html( *pre_stream);
	assert_file_write( *pre_stream, pre_main_filename);
	delete   pre_stream;
    } else
	cout << endl;

    assert_file_write( *index_stream, macroX( "\\lciIndexFilename"));
    delete index_stream;

    if (macroIsTrue("\\lciIfMultipleParts")) {
      *contents_stream << "<!-- End last part's contents table -->" << endl;
      *contents_stream << "</TABLE></TD></TABLE>" << endl; 
    }

    if ( ! noheader_switch)
	copy_and_filter_config_file( macroX( "\\lciTocFooter"),
				     *contents_stream);

    assert_file_write( *contents_stream, 
		       macroX( "\\lciContentsFilename"));
    delete contents_stream;

    assert_file_write( *anchor_stream, macroX( "\\lciAnchorFilename"));
    delete anchor_stream;

    return 0;
}

// EOF //

