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
#include <config.h>
#include <database.h>


/* Declarations from syntax.y */
/* ========================== */
extern int yydebug;
void yyparse();

/* Declarations from lex.yy   */
/* ========================== */
void init_scanner( FILE* in);


/* main */
/* ==== */
 
#define MaxParameters          1000
#define MaxOptionalParameters  1000
#define ErrParameters          10000
 
typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
Switch  trace_switch  = NO_SWITCH;
 
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
        cerr << "Usage: cgal_extact [-trace] [<infile> ...]" << endl;
        cerr << "       -trace       sets the `yydebug' variable of bison"
             << endl;
        exit(1);
    }

    if ( nParameters) {
        for ( i = 0; i < nParameters; i++) {
	    FILE* in;
            if ( (in = fopen( parameters[i], "r")) == NULL) {
	        fprintf( stderr, "\nbinsort: error: cannot open infile %s.\n",
                         parameters[0]);
                exit(1);
            }
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

