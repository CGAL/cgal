/**************************************************************************

  cc_extract_html.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : main program, command line parameter parsing
  System    : bison, flex, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <ctype.h>

#include <mstring.h>

#include <macro_dictionary.h>
#include <internal_macros.h>
#include <config.h>
#include <error.h>
#include <syntax.h>
#include <string_conversion.h>
#include <input.h>
#include <output.h>

using namespace std;


/* Program name and release          */
/* ================================= */
const string prog_name    = "cc_extract_html";
const string prog_release = "$Id$";

/* Configurable command line options */
/* ================================= */
Switch  V_switch           = NO_SWITCH;

Switch  trace_switch       = NO_SWITCH;
Switch  line_switch        = NO_SWITCH;

Switch  config_switch      = NO_SWITCH;
Switch  quiet_switch       = NO_SWITCH;
Switch  verbose_switch     = NO_SWITCH;
Switch  macro_def_switch   = NO_SWITCH;
Switch  macro_exp_switch   = NO_SWITCH;
Switch  macro_def2_switch  = NO_SWITCH;
Switch  macro_exp2_switch  = NO_SWITCH;
Switch  sty_macro_switch   = NO_SWITCH;
Switch  stack_trace_switch = NO_SWITCH;


/* Config filename:                  */
/* ================================= */
string sty_filename      = "latex.sty";


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
    insertInternalGlobalMacro( "\\lciExtractHtmlRelease", prog_release);
    insertInternalGlobalMacro( "\\lciManualDate",     manual_date);
    // check for date format as provided by latex_to_html
    if ( std::count( manual_date.begin(), manual_date.end(), ',') == 2) {
        // skip the day
        std::size_t pos = std::find( manual_date.begin(), manual_date.end(),
                                     ',') - manual_date.begin();
        if ( manual_date[pos+1] == ' ' )
            ++pos;
        string aux_date = manual_date;
        aux_date.replace( 0, pos, string());
        insertInternalGlobalMacro( "\\today", aux_date);

    } else {
        // else keep the full date as provided
        insertInternalGlobalMacro( "\\today", manual_date);
    }
    insertInternalGlobalMacro( "\\lciManualRelease",  manual_release);
    insertInternalGlobalMacro( "\\lciManualTitle",    manual_title);
    insertInternalGlobalMacro( "\\lciManualAuthor",   manual_author);
    insertInternalGlobalMacro( "\\lciIfVersionFlag",  V_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfQuietFlag",    quiet_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfVerboseFlag",  verbose_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfMacroDefFlag", macro_def2_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciIfMacroExpFlag", macro_exp2_switch
			       ? "\\lcTrue" : "\\lcFalse");
    insertInternalGlobalMacro( "\\lciOutputFilename", "<cout>");
    insertInternalGlobalMacro( "\\lciMainFilename",   "<cout>");
}


/* Index */
/* ======*/

extern int HREF_counter;



/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleString( const char* s) {
  if( current_ostream )
    *current_ostream << s;
}

void handleString( const string& s) {
  if( current_ostream )
    *current_ostream << s;
}

void handleChar( char c) {
  if( current_ostream )
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
    if( s )
       config_path = s;
    s = getenv("LATEX_CONV_INPUTS");
    if( s )
       latex_conv_inputs = s;

    int i;
    int nParameters = 0;
    typedef const char* ParamType;
    ParamType parameters[ MaxParameters + 1];

    Switch help_switch  = NO_SWITCH;
    Switch dummy_switch;

    insertInternalGlobalMacro( "\\lciInstallLatexConverterCSSFile", "" );

    for (i = 1; i < argc; i++) {
        /* check switches */
        detectSwitch( dummy_switch, "get_latex_conv_config");
            cout << config_path << endl;
            return 0;
        endDetect();
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
                if ( s[0] == '/') {
                    insertGlobalMacro( "\\lciIfRelativeHeaderPath",
                                       "<command line option>", 0,"\\lcFalse");
                } else {
                    insertGlobalMacro( "\\lciIfRelativeHeaderPath",
                                       "<command line option>", 0, "\\lcTrue");
                }
	    } else {
	        cerr << "*** Error: option -header needs an additional "
		        "parameter" << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "color");
            enableColor();
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
        detectSwitch( verbose_switch, "v");
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

        detectSwitch( V_switch, "V");
	    cerr << prog_name << " " << prog_release << " (c) Lutz Kettner"
		 << endl;
	    cerr << "Using: ";
            eraseMacro( "\\lciInstallLatexConverterCSSFile" );
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
        cerr << "       -quiet              no output." << endl;
        cerr << "       -v                  verbose, gives more context in "
                                           "error messages." << endl;
        cerr << "       -macrodef           trace macro definitions." << endl;
        cerr << "       -macroexp           trace macro expansions." << endl;
        cerr << "       -stymacro           trace style macros as well."<<endl;
        cerr << "       -stacktrace         stack trace for each error."<<endl;
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

    if (! include_stack.push_file( config_path + sty_filename))
	exit(1);
    yyparse();

    if (V_switch) {
	cerr << endl;
	exit(0);
    }


    main_stream   = &cout;

    // Prepare several streams:
    contents_stream = open_file_for_write( tmp_path +
					   macroX( "\\lciContentsFilename"));
    comments_stream = open_file_for_write( tmp_path + "comments.xml" );

    package_overview_stream = open_file_for_write( tmp_path +
                                              macroX( "\\lciPkgOverviewFilename") );

    index_stream = open_file_for_write( tmp_path +
					   macroX( "\\lciIndexFilename"));
    HREF_stream =  open_file_for_write( tmp_path +
                                           macroX( "\\lciHREFFilename"));

    HREF_counter =  open_counter_file_for_read( tmp_path +
                                           macroX("\\lciHREFCounterFilename"));

    if ( ! pre_main_filename.empty()) {
	pre_stream = open_file_for_write_with_path(tmp_path+pre_main_filename);
	current_filename  = pre_main_filename;
	pre_main_basename = basename_string( pre_main_filename);
	pre_main_rootname = rootname_string( pre_main_basename);
	current_basename  = pre_main_basename;
	current_rootname  = pre_main_rootname;
        pre_main_filepath = path_string( pre_main_filename);
        pre_main_uppath   = uppath_string( pre_main_filepath);
        current_filepath  = pre_main_filepath;
        current_uppath    = pre_main_uppath;
    }


    macro_def_switch = macro_def2_switch;
    macro_exp_switch = macro_exp2_switch;
    if ( trace_switch)
	yydebug = 1;


    for ( i = 0; i < nParameters; i++) {
	if ( ! pre_main_filename.empty()) {
	    main_stream      = pre_stream;
	    main_filename    = pre_main_filename;
	    main_basename    = pre_main_basename;
	    main_rootname    = pre_main_rootname;
	    main_filepath    = pre_main_filepath;
	    main_uppath      = pre_main_uppath;
	    insertInternalGlobalMacro( "\\lciMainFilename",  main_filename);
	    insertInternalGlobalMacro( "\\lciMainBasename",  main_basename);
	    insertInternalGlobalMacro( "\\lciMainRootname",  main_rootname);
	    insertInternalGlobalMacro( "\\lciMainPath",      main_filepath);
	    insertInternalGlobalMacro( "\\lciMainUppath",    main_filepath);
	} else {
	    main_stream      = &cout;
	}
        current_ostream  = main_stream;
        current_filename = main_filename;
        current_basename = main_basename;
        current_rootname = main_rootname;
        current_filepath = main_filepath;
        current_uppath   = main_uppath;
        insertInternalGlobalMacro( "\\lciOutputFilename",current_filename);
        insertInternalGlobalMacro( "\\lciOutputBasename",current_basename);
        insertInternalGlobalMacro( "\\lciOutputRootname",current_rootname);
        insertInternalGlobalMacro( "\\lciOutputPath",    current_filepath);
        insertInternalGlobalMacro( "\\lciOutputUppath",  current_uppath);


        anchor_stream = open_file_for_write( tmp_path +
                                             macroX( "\\lciAnchorFilename"));
        global_anchor_stream = anchor_stream;
        main_anchor_stream   = anchor_stream;

	if ( include_stack.push_file( parameters[i]))
	    yyparse();

	include_stack.push_string( "<end of conversion>",
				   "\\lciCheckNestingScopes",
				   0);
	include_stack.push_string( "<end of conversion>",
				   "\\lciEndOfConversion",
				   0);
	yyparse();


	assert_file_write( *main_stream, main_filename);

        if ( anchor_stream != 0 && anchor_stream != main_anchor_stream
             && anchor_stream != global_anchor_stream) {
            assert_file_write( *anchor_stream, macroX( "\\lciAnchorFilename"));
            delete anchor_stream;
            anchor_stream = 0;
        }
        if ( main_anchor_stream != 0
             && main_anchor_stream != global_anchor_stream) {
            assert_file_write( *main_anchor_stream,
                               macroX( "\\lciAnchorFilename"));
            delete main_anchor_stream;
            main_anchor_stream = 0;
        }
        if ( global_anchor_stream != 0) {
            assert_file_write( *global_anchor_stream,
                               macroX( "\\lciAnchorFilename"));
            delete global_anchor_stream;
            global_anchor_stream = 0;
        }

	if ( main_stream != &cout && main_stream != pre_stream) {
	    assert_file_write( *main_stream, main_filename);
	    delete   main_stream;
	    main_stream = &cout;
	    main_filename = "<cout>";
	    main_basename = main_filename;
	    main_rootname = main_filename;
	    main_filepath = string();
	    main_uppath = string();
	}
    }
    if ( ! quiet_switch)
	cerr << ']' << endl;

    if ( ! pre_main_filename.empty()) {
	assert_file_write( *pre_stream, pre_main_filename);
	delete   pre_stream;
    } else
	cout << endl;


    assert_file_write( *index_stream, macroX( "\\lciIndexFilename"));
    delete index_stream;

    assert_file_write( *package_overview_stream, macroX( "\\lciPkgOverviewFilename") );
    delete package_overview_stream;

    HREF_counter_stream = open_file_for_write( tmp_path +
                                 macroX( "\\lciHREFCounterFilename"));
    *HREF_counter_stream << HREF_counter;
    assert_file_write( *HREF_counter_stream,
                        macroX( "\\lciHREFCounterFilename"));
    delete HREF_counter_stream;


    assert_file_write( *contents_stream,
		       macroX( "\\lciContentsFilename"));
    delete contents_stream;


    assert_file_write( *HREF_stream, macroX( "\\lciHREFFilename"));
    delete HREF_stream;

    assert_file_write( *comments_stream, "comments.xml" );
    delete comments_stream;


    return firstError(); // reports non-zero return codes if there were
                         // non-aborting errors.
}

// EOF //

