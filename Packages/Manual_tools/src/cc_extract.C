/**************************************************************************
 
  cgal_extract.cc
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : main program, command line parameter parsing
  System    : bison, flex, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stream.h>
#include <database.h>
#include <config.h>

typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
Switch  trace_switch  = NO_SWITCH;
Switch  noc_switch    = NO_SWITCH;
Switch  nomc_switch   = NO_SWITCH;
Switch  nosc_switch   = NO_SWITCH;


/* Declarations from syntax.y */
/* ========================== */
extern int yydebug;
void yyparse();

/* Declarations from lex.yy   */
/* ========================== */
void init_scanner( FILE* in);


/* Name the scanned file */
/* ===================== */
const char* file_name = "<stdin>";


/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleMainComment( const Text& T) {
    if ( noc_switch || nomc_switch)
        return;
    if ( indentation_number() > 0) {
        cout << outdent;
	printComment( cout, T, true);
        cout << indent;
    } else {
	printComment( cout, T, true);
    }
}

void handleComment( const Text& T) {
    if ( noc_switch || nosc_switch)
        return;
    cout << indent;
    printComment( cout, T, false);
    cout << outdent;
}

void handleClass( const char* classname) {
    cout << newline;
    cout << "class " << classname << " {" << newline;
    cout << "public:" << newline;
    cout << indent;
}

void handleClassEnd( void) {
    cout << outdent;
    cout << newline;
    cout << "}" << newline;
}

void handleClassTemplate( const char* classname) {
    cout << newline;
    cout << "template < class ";
    const char* s = classname;
    while ( *s != 0 && *s != '<') s++;
    if ( *s == 0)
        printErrorMessage( TemplateParamExpectedError);
    else {
        int nesting = 0;
	s++;
	while ( nesting >= 0 && *s != 0) {
	    switch ( *s) {
	    case '<':
	        nesting ++;
		cout << *s;
		break;
	    case '>':
	        nesting --;
		if ( nesting >= 0)
		    cout << *s;
		break;
	    case ',':
	        if ( nesting == 0)
		    cout << ", class ";
		else
		    cout << *s;
		break;
	    default:
		cout << *s;
		break;
	    }
	    s++;
	}
	if ( nesting >= 0)
	    printErrorMessage( MalformedTemplateParamError);
    }
    cout << " >" << newline;
    cout << "class ";
    s = classname;
    while ( *s != 0 && *s != '<') {
        cout << *s;
	s++;
    }
    cout << " {" << newline;
    cout << "public:" << newline;
    cout << indent;
}

void handleClassTemplateEnd( void) {
    cout << outdent;
    cout << newline;
    cout << "}" << newline;
}


void handleDeclaration( const char* decl) {
    cout << endl << newline;
    cout << decl;
}

void handleFunctionDeclaration( const char* decl) {
    cout << endl << newline;
    cout << decl;
}

void handleFunctionTemplateDeclaration( const char* templ, const char* decl) {
    cout << endl << newline;
    cout << "template < class ";
    const char* s = templ;
    int nesting = 0;
    while ( *s != 0) {
	    switch ( *s) {
	    case '<':
	        nesting ++;
		cout << *s;
		break;
	    case '>':
	        nesting --;
		if ( nesting >= 0)
		    cout << *s;
		break;
	    case ',':
	        if ( nesting == 0)
		    cout << ", class ";
		else
		    cout << *s;
		break;
	    default:
		cout << *s;
		break;
	    }
	    s++;
    }
    if ( nesting != 0)
        printErrorMessage( MalformedTemplateParamError);
    cout << " >" << newline;
    cout << decl;
}


/* main */
/* ==== */
 
#define MaxParameters          1000
#define MaxOptionalParameters  1000
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
    int i;
    int nParameters = 0;
    char *parameters[ MaxParameters + 1];
 
    Switch help_switch = NO_SWITCH;
 
    for (i = 1; i < argc; i++) {

        /* check switches */
        detectSwitch( trace_switch, "trace");
	    yydebug = 1;
        endDetect();

        detectSwitch( noc_switch, "noc");
        endDetect();
        detectSwitch( nomc_switch, "nomc");
        endDetect();
        detectSwitch( nosc_switch, "nosc");
        endDetect();

        detectSwitch( help_switch, "h");
        endDetect();
        detectSwitch( help_switch, "H");
        endDetect();
        detectSwitch( help_switch, "help");
        endDetect();
 
        /* else get standard or optional paramters */
        if ( nParameters < MaxParameters ) {
            parameters[nParameters ++] = argv[i];
            continue;
        }
 
        nParameters = ErrParameters;
        break;
    }
 
    if ((nParameters < MaxParameters - MaxOptionalParameters) ||
        (nParameters > MaxParameters) || (help_switch != NO_SWITCH)) {
        if (help_switch == NO_SWITCH)
            cerr << "Error: in parameter list" << endl;
        cerr << "Usage: cgal_extract [<options>] [<infile> ...]" << endl;
        cerr << "       -nomc        no main comments"  << endl;
        cerr << "       -nosc        no sub comments"  << endl;
        cerr << "       -noc         no comments (main and sub)"  << endl;
        cerr << "       -trace       sets the `yydebug' variable of bison"
             << endl;
        exit(1);
    }

    if ( nParameters) {
        for ( i = 0; i < nParameters; i++) {
	    FILE* in;
            if ( (in = fopen( parameters[i], "r")) == NULL) {
	        fprintf( stderr, 
			 "\ncgal_extract: error: cannot open infile %s.\n",
                         parameters[i]);
                exit(1);
            }
	    file_name = parameters[i];
	    init_scanner( in);
	    yyparse();
            fclose( in);
	}
    } else {
        init_scanner( stdin);
	yyparse();
    }
    cout << endl;
    return 0;
}

// EOF //

