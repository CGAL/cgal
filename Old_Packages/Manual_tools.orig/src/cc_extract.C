/**************************************************************************
 
  cc_extract.C
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Extract C++ declarations from a TeX specification
              file written with the cc_manual.sty.
              Main program, command line parameter parsing.
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
#include <buffer.h>
#include <config.h>

// Max. width of the formatted output
#define MaxTextWidth 74

typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
Switch  trace_switch  = NO_SWITCH;
Switch  noc_switch    = NO_SWITCH;
Switch  nomc_switch   = NO_SWITCH;
Switch  nosc_switch   = NO_SWITCH;

// Format comments
// ==========================
int printComment( ostream &out, const Buffer_list& T, bool leadingLine) {
    int  width = MaxTextWidth - indentation_number() - 3;
    int  w     = width;
    int  state = 0;     // 0 = start, 1 = after token, 2 = spaces
                        // 3 = one newline (and spaces), 4 = newlines
    for ( Buffer_const_iterator words = T.begin(); words != T.end(); ++words) {
	bool is_space = (*words)->is_space();
	int  len      = int((*words)->size()) - 1;
	switch ( state) {
	case 0:
	    if ( ! is_space) {
		if ( leadingLine)
		    out << endl;
		out << ind_newline << "// " << (*words)->begin();
		w -= len;
		state = 1;
	    }
	    break;
	case 1:
	    if ( ! is_space || len > 0) {
		if ( ! is_space) {
		    if ((len > w) && (w != width)) {
			w = width;
			out << ind_newline << "// ";
		    }
		    out << (*words)->begin();
		    w -= len;
		} else {
		    if ( (*words)->begin()[0] == '\n')
			state = 3;
		    else
			state = 2;
		}
	    }
	    break;
	case 2:
	    if ( ! is_space || len > 0) {
		if ( ! is_space) {
		    if ((len >= w) && (w != width)) {
			w = width;
			out << ind_newline << "// ";
		    } else {
			out << ' ';
			w--;
		    }
		    out << (*words)->begin();
		    w -= len;
		    state = 1;
		} else {
		    if ( (*words)->begin()[0] == '\n')
			state = 3;
		}
	    }
	    break;
	case 3:
	    if ( ! is_space || len > 0) {
		if ( !is_space) {
		    if ((len >= w) && (w != width)) {
			w = width;
			out << ind_newline << "// ";
		    } else {
			out << ' ';
			w--;
		    }
		    out << (*words)->begin();
		    w -= len;
		    state = 1;
		} else {
		    if ( (*words)->begin()[0] == '\n')
			state = 4;
		}
	    }
	    break;
	case 4:
	    if ( ! is_space || len > 0) {
		if ( ! is_space) {
		    out << ind_newline << "// ";
		    out << ind_newline << "// ";
		    w = width;
		    out << (*words)->begin();
		    w -= len;
		    state = 1;
		}
	    }
	    break;
	}
    }
    return state;
}

/* Declarations from syntax.y */
/* ========================== */
extern int yydebug;
void yyparse();

/* Declarations from lex.yy   */
/* ========================== */
void init_scanner( FILE* in);
extern char*   global_classname;
extern char*   global_template_params;
extern char*   global_ref_name;


/* Name the scanned file */
/* ===================== */
const char* file_name = "<stdin>";


/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleMainComment( const Buffer_list& T) {
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

void handleComment( const Buffer_list& T) {
    if ( noc_switch || nosc_switch)
        return;
    cout << indent;
    printComment( cout, T, false);
    cout << outdent;
}

void handleClass( const char* classname) {
    if ( global_classname)
        free( global_classname);
    global_classname = strdup( classname);
    const char* s = classname;
    while ( *s != 0 && *s != '<') 
	s++;

    if ( *s == 0) {
	cout << ind_newline;
	cout << "class " << classname;
    } else {
	global_template_params = global_classname;
	while( *global_template_params && *global_template_params != '<')
	    ++global_template_params;
	cout << ind_newline;
	cout << "template < class ";
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
	cout << " >" << ind_newline;
	cout << "class ";
	s = classname;
	while ( *s != 0 && *s != '<') {
	    cout << *s;
	    s++;
	}
    }
    cout << " {" << ind_newline;
    cout << "public:" << ind_newline;
    cout << indent;
}

void handleClassEnd( void) {
    global_template_params = 0;
    if ( global_classname)
        free( global_classname);
    global_classname = NULL;
    cout << outdent;
    cout << ind_newline;
    cout << "};" << ind_newline;
}

void handleRefPage( const char* token) {
    if ( global_ref_name)
        free( global_ref_name);
    global_ref_name = strdup( token);
}

void handleRefPageEnd( void) {
    if ( global_ref_name)
        free( global_ref_name);
    global_ref_name = NULL;
}

void handleDeclaration( const char* decl) {
    cout << endl << ind_newline;
    cout << decl;
}

void handleNestedType( const char* decl) {
    cout << endl << ind_newline;
    cout << "Nested type required: " << global_classname << "::" << decl;
}

void handleMethodDeclaration( const char* decl) {
    cout << endl << ind_newline;
    cout << decl;
}

void handleFunctionDeclaration( const char* decl) {
    cout << endl << ind_newline;
    cout << decl;
}

void handleFunctionTemplateDeclaration( const char* templ, const char* decl) {
    cout << endl << ind_newline;
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
    cout << " >" << ind_newline;
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
 
int main( int argc, char **argv) {
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
	cerr << "cc_extract $Revision$ (c) Lutz Kettner" << endl;
        cerr << "Usage: cc_extract [<options>] [<infile> ...]" << endl;
        cerr << "       -nomc        no main comments"  << endl;
        cerr << "       -nosc        no sub comments"  << endl;
        cerr << "       -noc         no comments (main and sub)"  << endl;
        cerr << "       -trace       sets the `yydebug' variable of bison"
             << endl;
        cerr << "       -h, -help    this help message"  << endl;
        exit(1);
    }

    if ( nParameters) {
        for ( i = 0; i < nParameters; i++) {
	    FILE* in;
            if ( (in = fopen( parameters[i], "r")) == NULL) {
	        fprintf( stderr, 
			 "\ncc_extract: error: cannot open infile %s.\n",
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

