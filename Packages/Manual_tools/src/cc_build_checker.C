/**************************************************************************
 
  cc_build_checker.cc
  =============================================================
  Project   : CGAL tool that constructs the specification checker
  Function  : Checks whether the C++ declarations from a TeX 
              specification file written with the cc_manual.sty
              really exists in a C++ header file. This program 
              extracts a flex program from a TeX specification
              that does the actual checking. This program is used
              with a script `cc_check' that is provided for the user.
              main program, command line parameter parsing
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
#include <ctype.h>
#include <fstream.h>
#include <buffer.h>
#include <config.h>

#define MaxParameters          1002
#define MaxOptionalParameters  1000
#define ErrParameters          10000

int nParameters = 0;
char *parameters[ MaxParameters + 1];

typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
Switch  trace_switch  = NO_SWITCH;


/* Declarations from syntax.y */
/* ========================== */
extern int yydebug;
void yyparse();

/* Declarations from lex.yy   */
/* ========================== */
void init_scanner( FILE* in);
extern int line_number;
extern int unchecked_tag;
extern char*   global_classname;
extern char*   global_template_params;
extern char*   global_ref_name;

/* Name the scanned file */
/* ===================== */
const char* file_name = "<stdin>";
const char* error_filename;

/* Global variables                            */
/* =========================================== */
/* Buffer to hold an appropriate error message */
Buffer* errmessage = new Buffer;

/* the n_class that counts the actual class */
int n_class = 0;


/* Buffer to hold the pattern to be scanned for */
/* Only necessary for special patterns */
Buffer* pattern  = new Buffer;
Buffer* pattern2 = new Buffer;

/* Buffer to hold a generalized function declaration where 
   in the parameter lists variable names are substituted by the
   special constant IDFIER_TAG
   */
Buffer* function_decl  = new Buffer;
const char IDFIER_TAG = 254;
const char INLINE_TAG = 253;

/* Number of pattern's written */
int n_pattern = 0;

ostream* p_out_pattern;
ostream* p_out_errmessage;


/* Functions to write the checker program */
/* ====================================== */
char char_tmp[2]         = " ";
char char_tricky_tmp[16] = "\"{ws}*\" \"{ws}*\"";

// const char* patternText( char c, char d = 'X') {
const char* patternText( char c, char) {
    switch ( c) {
    case IDFIER_TAG:
        return "\"{idfier}\"";
    case INLINE_TAG:
        return "\"{ws}*\">\"{ws}*(inline)?{ws}*\"";
    case ' ':
    case '\t':
    case '\r':
    case '\n':
        return "\"{ws}*\"";
    case '\\':
        return "\\\\";
    case '"':
        return "\\\"";
    case ',':    // A couple of cases where whitespaces can follow
    case '{':
    case '}':
    case '=':
    case '(':
    case ')':
    case '<':
    case '>':
    case '*':
    case '&':
          // if ( d != ' ') {
	char_tricky_tmp[7] = c;
        return char_tricky_tmp;
	  // }
	// break;
    }
    char_tmp[0] = c;
    return char_tmp;
}

void printProlog( ostream& out_pattern, ostream& out_errmessage) {
    int i;
    out_pattern << "/* =============================================" << endl;
    out_pattern << "   CC checker program to verify consistency" << endl;
    out_pattern << "   between specification and implementation." << endl;
    out_pattern << "   Generated for the specification:" << endl;
    out_pattern << "       `" << file_name;
    for ( i=3; i<nParameters; i++) {
        out_pattern << " " << parameters[ i];
    }
    out_pattern << "'" << endl;
    out_pattern << endl;
    out_pattern << "   R 1.2    03.12.1996   Lutz Kettner" << endl;
    out_pattern << "            kettner@acm.org" << endl;
    out_pattern << "            ETH Zurich, Switzerland" << endl;
    out_pattern << endl;
    out_pattern << "   Compilation with `flex' and `cc' (or `gcc')" << endl;
    out_pattern << "============================================== */" << endl;
    out_pattern << endl;
    out_pattern << "%{" << endl;
    out_pattern << "#include <stdlib.h>" << endl;
    out_pattern << "#include <stdio.h>" << endl;
    out_pattern << "#include <string.h>" << endl;
    out_pattern << "#include \"" << error_filename << "\"" << endl;
    out_pattern << endl;
    out_pattern << "/* the tag array that collects correct matches */" << endl;
    out_pattern << "int *tags;" << endl;
    out_pattern << endl;
    out_pattern << "/* the n_class that counts the actual class */" << endl;
    out_pattern << "int n_class = -1;" << endl;
    out_pattern << endl;
    out_pattern << "/* the name array that stores the filename */" << endl;
    out_pattern << "const char **name;" << endl;
    out_pattern << endl;
    out_pattern << "/* mention filename */" << endl;
    out_pattern << "const char* file_name = \"<stdin>\";" << endl;
    out_pattern << endl;
    out_pattern << "/* count linenumbers */" << endl;
    out_pattern << "int line_number = 1;" << endl;
    out_pattern << endl;
    out_pattern << "/* command line options */" << endl;
    out_pattern << "int funnel   = 0;" << endl;
    out_pattern << "int match    = 0;" << endl;
    out_pattern << "int unmatch  = 0;" << endl;
    out_pattern << "int comm     = 0;" << endl;
    out_pattern << endl;
    out_pattern << "/* Make the debug options easier */" << endl;
    out_pattern << "#define MD if( match) fputs( yytext, stdout)" << endl;
    out_pattern << "#define UD if( unmatch) fputs( yytext, stdout)" << endl;
    out_pattern << "#define CD if( comm) fputs( yytext, stdout)" << endl;
    out_pattern << "#define AD if( match || unmatch || comm) "
                   "fputs( yytext, stdout)" << endl;
    out_pattern << endl;
    out_pattern << "/* keep track of newlines in matched lines */" << endl;
    out_pattern << "void count_newlines( char* s){" << endl;
    out_pattern << "    while (*s != 0) {" << endl;
    out_pattern << "        if ( *s == '\\n')" << endl;
    out_pattern << "            line_number++;" << endl;
    out_pattern << "        s++;" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "}" << endl;
    out_pattern << endl;
    out_pattern << "/* Hack, to get rid of the yywrap. */" << endl;
    out_pattern << "#define YY_SKIP_YYWRAP 1" << endl;
    out_pattern << "#define yywrap() 1" << endl;
    out_pattern << "" << endl;
    out_pattern << "%}" << endl;
    out_pattern << endl;
    out_pattern << "%x COMMENT" << endl;
    out_pattern << "%x STRING" << endl;
    out_pattern << "%x FUNNEL" << endl;
    out_pattern << endl;
    out_pattern << "/* A couple of abbreviations */" << endl;
    out_pattern << "cccomment            \"//\".*$" << endl;
    out_pattern << "cccommentnl          \"//\".*[\\n]" << endl;
    out_pattern << "fwcomment            [@][!].*$" << endl;
    out_pattern << "ws                   [ \\t\\n\\r]|{cccommentnl}" << endl;
    out_pattern << "nol                  [^a-zA-Z0-9_]" << endl;
    out_pattern << "classend             {ws}*(([:][^:])|([{]))" << endl;
    out_pattern << "optionalclass        {ws}*(class{ws}*)?" << endl;
    out_pattern << "optionalkeywords     ([a-z]+({ws}+[a-z]+)*)?" << endl;
    out_pattern << "idfier               [a-zA-Z_][a-zA-Z0-9_]*" << endl;
    out_pattern << endl;
    out_pattern << "%%" << endl;
    out_pattern << "    /* The matching rules */" << endl;
    out_pattern << "    /* ================== */" << endl;
    out_pattern << "    if ( funnel) {" << endl;
    out_pattern << "        BEGIN( FUNNEL);" << endl;
    out_pattern << "    }" << endl;
    out_pattern << endl;
    out_pattern << "<INITIAL,COMMENT,STRING,FUNNEL>[\\n]   "
                   "{ line_number ++; AD;}" << endl;
    out_pattern << endl;
    out_pattern << "{cccomment}               { CD; /* ignore it */ }" << endl;
    out_pattern << "\"/*\"                      { CD; BEGIN( COMMENT); }" 
		<< endl;
    out_pattern << "<COMMENT>\"*/\"             { CD; BEGIN( INITIAL); }" 
		<< endl;
    out_pattern << "<COMMENT>.                { CD; /* ignore it */ }" << endl;
    out_pattern << endl;
    out_pattern << "\"'\\\\''\"                   { CD;/* ignore */ }" << endl;
    out_pattern << "\"'\".*\"'\"                  { CD;/* ignore */ }" << endl;
    out_pattern << "<INITIAL,STRING>\"\\\\\\\"\"    { CD; /* ignore */ }"
		<< endl;
    out_pattern << "\"\\\"\"                      { CD; BEGIN( STRING); }" 
		<< endl;
    out_pattern << "<STRING>\"\\\"\"              { CD; BEGIN( INITIAL);}" 
		<< endl;
    out_pattern << "<STRING>.                 { CD; /* ignore */ }" << endl;
    out_pattern << endl;
    out_pattern << "<FUNNEL>\"@{\"|\"@begin\"    { CD; BEGIN( INITIAL); }" 
		<< endl;
    out_pattern << "\"@}\"|\"@end\"              { CD; if ( funnel) "
                   "BEGIN( FUNNEL); }" << endl;
    out_pattern << "<FUNNEL>{fwcomment}       { CD; /* ignore */ }" << endl;
    out_pattern << "<FUNNEL>.                 { CD; /* ignore */ }" << endl;
    out_pattern << endl;

    out_errmessage << "/* =============================================" 
		   << endl;
    out_errmessage << "   CC checker program to verify consistency" << endl;
    out_errmessage << "   between specification and implementation." << endl;
    out_errmessage << "   Generated for the specification:" << endl;
    out_errmessage << "       `" << file_name;
    for ( i=3; i<nParameters; i++) {
        out_errmessage << " " << parameters[ i];
    }
    out_errmessage << "'" << endl;
    out_errmessage << endl;
    out_errmessage << "   R 1.2    03.12.1996   Lutz Kettner" << endl;
    out_errmessage << "            kettner@acm.org" << endl;
    out_errmessage << "            ETH Zurich, Switzerland" << endl;
    out_errmessage << endl;
    out_errmessage << "   Auxiliary file for the pattern texts." << endl;
    out_errmessage << "============================================= */" 
		   << endl;
    out_errmessage << endl;
    out_errmessage << "const char* patternText( int number) {" << endl;
    out_errmessage << "    switch ( number) {" << endl;
}

void printEpilog( ostream& out_pattern, ostream& out_errmessage) {
    out_pattern << ".               {  UD; /* ignore inbetween */ }" << endl;
    out_pattern << endl;
    out_pattern << "%%" << endl;
    out_pattern << endl;
    out_pattern << "int main( int argc, const char* argv[]) {" << endl;
    out_pattern << "    int i;" << endl;
    out_pattern << "    int err = 0;" << endl;
    out_pattern << "    int flags = 1;" << endl;
    out_pattern << "    int err_count = 0;" << endl;
    out_pattern << "    int warn_count = 0;" << endl;
    out_pattern << "    tags = (int*)malloc( " << n_pattern 
		<< " * sizeof( int));" << endl;
    out_pattern << "    name = (const char**)malloc( " << n_pattern 
		<< " * sizeof( char*));" << endl;
    out_pattern << "    if ( ! tags || ! name) {" << endl;
    out_pattern << "        fprintf( stderr, \"fatal error: " 
		   "not enough memory to allocate " << n_pattern 
		<< " int's and char**'s.\\n\"); "<< endl;
    out_pattern << "        exit( 1);" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    memset( tags, 0, " << n_pattern << " * sizeof( int));"
		<< endl;
    out_pattern << "    if ( argc > 1 && strcmp( argv[1], \"-h\") == 0) {" 
		<< endl;
    out_pattern << "        fprintf( stderr, \"Usage: %s "
                   "[<options>] [<impl-files> ...]\\n\", argv[0]);" << endl;
    out_pattern << "        fprintf( stderr, \"           -h         "
                   "help\\n\");" << endl;
    out_pattern << "        fprintf( stderr, \"           -fw        "
                   "recognize FunnelWeb or AnyWeb keywords\\n\");" << endl;
    out_pattern << "        fprintf( stderr, \"           -match     "
                   "print all matched text\\n\");" << endl;
    out_pattern << "        fprintf( stderr, \"           -unmatch   "
                   "print all unmatched text\\n\");" << endl;
    out_pattern << "        fprintf( stderr, \"           -comm      "
                   "print all commented text\\n\");" << endl;
    out_pattern << "        exit( 1);" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( argc > flags && strcmp( argv[flags], \"-fw"
                   "\") == 0) {" << endl;
    out_pattern << "        funnel = 1;" << endl;
    out_pattern << "        flags++;" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( argc > flags && strcmp( argv[flags], \"-match"
                   "\") == 0) {" << endl;
    out_pattern << "        match = 1;" << endl;
    out_pattern << "        flags++;" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( argc > flags && strcmp( argv[flags], \"-unmatch"
                   "\") == 0) {" << endl;
    out_pattern << "        unmatch = 1;" << endl;
    out_pattern << "        flags++;" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( argc > flags && strcmp( argv[flags], \"-comm"
                   "\") == 0) {" << endl;
    out_pattern << "        comm = 1;" << endl;
    out_pattern << "        flags++;" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( argc > flags) {" << endl;
    out_pattern << "        for ( i=flags; i<argc; i++) {" << endl;
    out_pattern << "            FILE *in;" << endl;
    out_pattern << "            if ( (in = fopen( argv[i], \"r\")) == NULL) {"
		<< endl;
    out_pattern << "                fprintf( stderr," << endl;
    out_pattern << "                         \"error: cannot open file `%s' "
                   "for reading.\\n\"," << endl;
    out_pattern << "                         argv[i]);" << endl;
    out_pattern << "                exit( 1);" << endl;
    out_pattern << "            }" << endl;
    out_pattern << "            file_name = argv[i];" << endl;
    out_pattern << "            line_number = 1;" << endl;
    out_pattern << "            yyrestart( in);" << endl;
    out_pattern << "            yylex();" << endl;
    out_pattern << "            fclose( in);" << endl;
    out_pattern << "        }" << endl;
    out_pattern << "    } else {" << endl;
    out_pattern << "        line_number = 1;" << endl;
    out_pattern << "        yyrestart( stdin);" << endl;
    out_pattern << "        yylex();" << endl;
    out_pattern << "    }" << endl;
    out_pattern << endl;
    out_pattern << "    for ( i=0; i<" << n_pattern << "; i++) {" << endl;
    out_pattern << "        if ( tags[i] == 0) {" << endl;
    out_pattern << "            err = 2;" << endl;
    out_pattern << "            err_count++;" << endl;
    out_pattern << "            fprintf( stderr, \"check error: %s cannot "
                   "be found in the implementation.\\n\", patternText( i));" 
		<< endl;
    out_pattern << "        } else if ( tags[i] > 0) {" << endl;
    out_pattern << "            err = 3;" << endl;
    out_pattern << "            warn_count++;" << endl;
    out_pattern << "            fprintf( stderr, \"check warning: %s only "
                   "recognized as prefix in `%s' line %d.\\n\", "
                   "patternText( i), name[i], tags[i]);" << endl;
    out_pattern << "        }" << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    if ( err_count || warn_count) {" << endl;
    out_pattern << "        fprintf( stderr, \"---------------------"
                   "----------------------------------\\n\");  " << endl;
    out_pattern << "        fprintf( stderr, \"%d error(s) and %d warning(s)"
                   ".\\n\\n\", err_count, warn_count); " << endl;
    out_pattern << "        " << endl;
    out_pattern << "    }" << endl;
    out_pattern << "    free( tags);" << endl;
    out_pattern << "    free( name);" << endl;
    out_pattern << "    return err;" << endl;
    out_pattern << "}" << endl;
    out_pattern << "/* EOF */" << endl;

    out_errmessage << "    }" << endl;
    out_errmessage << "    fprintf( stderr, \"fatal internal error: unknown "
                      "pattern number %d occured.\", number);" << endl;
    out_errmessage << "    exit( 1);" << endl;
    out_errmessage << "    return \"\";" << endl;
    out_errmessage << "}" << endl;
    out_errmessage << "/* EOF */" << endl;
}

// Simplify the pattern: remove redundent sequences of "" and white spaces
//------------------------------------------------------------------------
void printReducedPattern( ostream& out_pattern, const char* s, int l) {
    // checks for multiple {ws}* and "" occurences.
    bool ws_plus = false;
    if ( l > 2 && s[l-1] == '"' && s[l-2] == '"')
        l -= 2;
    for ( int i=0; i<l; i++) {
        if ( i<l-1 && s[i] == '"' && s[i+1] == '"')
	    i ++;
	else if ( i+11 < l
		  && s[i] == '{'
		  && s[i+1] == 'w'
		  && s[i+2] == 's'
		  && s[i+3] == '}'
		  && s[i+5] == '"'
		  && s[i+6] == '"'
		  && s[i+7] == '{'
		  && s[i+8] == 'w'
		  && s[i+9] == 's'
		  && s[i+10] == '}'
		  && ( s[i+4]  == '*' || s[i+4]  == '+')
		  && ( s[i+11] == '*' || s[i+11] == '+'))
	{
	    if ( s[i+4] == '+' && ! ws_plus) {
	        out_pattern << "{ws}+";
		ws_plus = true;
	    }
	    i += 6;
	} else if ( i+4 < l
		  && s[i] == '{'
		  && s[i+1] == 'w'
		  && s[i+2] == 's'
		  && s[i+3] == '}'
		  && ( s[i+4]  == '*' || s[i+4]  == '+'))
	 {
	     if ( i+5 < l && ! ws_plus) {
	         out_pattern << "{ws}" << s[i+4];
		 ws_plus = (s[i+4] == '+');
	     }
	     i += 4;
	 } else {
	    out_pattern << s[i];
	    ws_plus = false;
	}
    }
}

// Print a simple declaration, results in a single matching rule
// -------------------------------------------------------------
void printPattern( ostream& out_pattern, 
		   ostream& out_errmessage,
		   Buffer *err,
		   Buffer *pat = 0) {
    if ( unchecked_tag) {
        unchecked_tag = 0;
	return;
    }
    int i;
    int l = err->size() - 1;
    const char *str = err->begin();
    if ( pat == 0) {
        pat = pattern;
	pat->flush();
	pat->add( '"');
	for ( i=0; i<l; i++)
	    pat->add( patternText( str[i], str[i+1]));
	pat->add( '"');
    }
    printReducedPattern( out_pattern, pat->begin(), pat->size() - 1);
    out_pattern << "  {" << endl;
    out_pattern << "        tags[" << n_pattern << "] = -1;" << endl;
    out_pattern << "        count_newlines( yytext);" << endl;
    out_pattern << "        MD;" << endl;
    out_pattern << "    }" << endl;
    out_errmessage << "    case " << n_pattern << ":" << endl;
    out_errmessage << "        return \"line " << line_number << " in `"
		   << file_name << "': \\\"";
    for ( i=0; i<l; i++) {
        if ( str[i] == '\\')
	    out_errmessage << "\\\\";
	else if ( str[i] == '"')
	    out_errmessage << "\\\"";
	else 
	    out_errmessage << str[i];
    }
    out_errmessage << "\\\"\";" << endl;
    n_pattern++;
} 

// Print a function declaration, results in two matching rules,
// one for the exact match and one for possibly trailing
// default parameters.
// -------------------------------------------------------------
void printFunctionPattern( ostream& out_pattern, 
			   ostream& out_errmessage,
			   Buffer *err,
			   Buffer *pat  = 0,
			   Buffer *pat2 = 0) {
    if ( unchecked_tag) {
        unchecked_tag = 0;
	return;
    }
    int  i;
    bool default_params = false;
    bool nonzero_param_list = false;
    int  nesting = 0;
    int  l = err->size() - 1;
    const char *str = err->begin();
    if ( pat == 0) {
        pat = pattern;
	pat->flush();
	pat->add( '"');
	for ( i=0; i<l; i++) {
	    switch ( str[ i]) {
	    case '(':
	        nesting ++;
		break;
	    case ')':
	        nesting --;
		if ( nesting == 0) {
		    default_params = true;
		    pat->add( "\"{ws}*\"");
		    pat2 = pattern2;
		    pat2->flush();
		    pat2->add( * pat);
		    if ( nonzero_param_list)
		        pat2->add( ',');
		    pat2->add( "\"{ws}*[^ \\n\\r\\t)]");
		}
		break;
	    default:
	        if ( nesting > 0 && str[i] > ' ')
		    nonzero_param_list = true;
	    }
	    pat->add( patternText( str[i], str[i+1]));
	    //	    if ( default_params)
	    //        pat2->add( patternText( str[i], str[i+1]));
	}
	pat->add( '"');
	// if ( default_params)
	//     pat2->add( '"');
    }
    printReducedPattern( out_pattern, pat->begin(), pat->size() - 1);
    out_pattern << "    {" << endl;
    out_pattern << "        tags[" << n_pattern << "] = -1;" << endl;
    out_pattern << "        count_newlines( yytext);" << endl;
    out_pattern << "        MD;" << endl;
    out_pattern << "    }" << endl;
    if ( default_params || pat2 != 0) {
        printReducedPattern( out_pattern, pat2->begin(), pat2->size() - 1);
        out_pattern << "    {" << endl;
	out_pattern << "        if (tags[" << n_pattern << "] == 0) {" << endl;
	out_pattern << "            tags[" << n_pattern << "] = line_number;"
		    << endl;
	out_pattern << "            name[" << n_pattern << "] = file_name;"
		    << endl;
	out_pattern << "        }" << endl;
	out_pattern << "        count_newlines( yytext);" << endl;
	out_pattern << "        MD;" << endl;
	out_pattern << "    }" << endl;
    }
    out_errmessage << "    case " << n_pattern << ":" << endl;
    out_errmessage << "        return \"line " << line_number << " in `"
		   << file_name << "': \\\"";
    for ( i=0; i<l; i++) {
        switch (str[i]) {
	case '\\':
	    out_errmessage << "\\\\";
	    break;
	case '"':
	    out_errmessage << "\\\"";
	    break;
	case IDFIER_TAG:
	    out_errmessage << "{idfier}";
	    break;
	case INLINE_TAG:
	    out_errmessage << ">";
	    break;
	default:
	    out_errmessage << str[i];
	    break;
	}
    }
    out_errmessage << "\\\"\";" << endl;

    n_pattern++;
} 



/* Translate argument names in function declarations to idfiers */
/* ============================================================ */

// Given a string interval, decide whether it contains a valid type name
// or only declarative keywords and other symbols. This check is necessary
// to distinguish between function arguments with and without argument
// name.
bool isTypeInFunctionArgument( const char* s, const char* end) {
    while( s != end) {
        // isolate identifier
        while( s != end && ( !isalpha(*s) && *s != '_'))
            s++;
	const char* p = s;  // idfier start
        while( s != end && ( isalpha(*s) || *s == '_'))
            s++;
	if ( s == p)
	    return false;
	bool is_keyword = false;
	is_keyword = is_keyword || 0 ==  strncmp( p, "auto", 4);
	is_keyword = is_keyword || 0 ==  strncmp( p, "class", 5);
	is_keyword = is_keyword || 0 ==  strncmp( p, "const", 5);
	is_keyword = is_keyword || 0 ==  strncmp( p, "enum", 4);
	is_keyword = is_keyword || 0 ==  strncmp( p, "extern", 6);
	is_keyword = is_keyword || 0 ==  strncmp( p, "register", 8);
	is_keyword = is_keyword || 0 ==  strncmp( p, "signed", 6);
	is_keyword = is_keyword || 0 ==  strncmp( p, "static", 6);
	is_keyword = is_keyword || 0 ==  strncmp( p, "struct", 6);
	is_keyword = is_keyword || 0 ==  strncmp( p, "template", 8);
	is_keyword = is_keyword || 0 ==  strncmp( p, "union", 5);
	is_keyword = is_keyword || 0 ==  strncmp( p, "unsigned", 8);
	is_keyword = is_keyword || 0 ==  strncmp( p, "volatile", 8);
	if (!is_keyword)
	    return true;
    }
    return false;
}

// Scan a full C++ function declaration in decl. Return value in func.
// Scans argument by argument and replaces the argument name by a symbol 
// denoting an idfier. The output routine replaces this symbol with the 
// actual regular expression for an identifier.
void translateFunctionArgumentLists( const char* decl, Buffer *func) {
    func->flush();
    const char* current = decl;
    while ( *current != '(') {
        current++;
        CC_Assert( *current);
    }
    CC_Assert( *current == '(');
    current++;
    CC_Assert( *current);
    func->add( decl, current - decl);

    /* the parsing starts here. It counts parenthesis. */
    int nesting = 1;
    bool in_init = false;
    const char* from = current;
    const char* last_idfier = current - 1;
    while ( nesting > 0) {
        CC_Assert( *current);
	if ( nesting == 1 && *current == '=')
	    in_init = true;
	if (      nesting == 1 
	       && !in_init 
	       && (isalnum(*current) || *current == '_'))
	    last_idfier = current;

	if ( nesting == 1 && (*current == ',' || *current == ')')) {
            if ( last_idfier - from < 0) { 
                /* fail save behavior in case of errors */
	        func->add( from, current-from+1);
	    } else {
		const char* p = last_idfier;
		while( isalnum( *p) || *p == '_')
		    p--;
		p++;  /* points to idfier start */
		CC_Assert( p - from >= 0);
		if ( isTypeInFunctionArgument( from, p)) {
		    func->add( from, p - from);
		    func->add( IDFIER_TAG);
		    func->add( last_idfier + 1, current - last_idfier);
		} else {
		    /* fail save behavior in case of errors */
		    func->add( from, current-from+1);
		}
	    }
	    from = current + 1;
	    if (*current == ')')
	        nesting--;
	    in_init = false;
	} else {
	    switch ( *current) {
	    case '(':
	    case '<':
	        nesting++;
		break;
	    case ')':
	    case '>':
	        nesting--;
		CC_Assert( nesting >= 1);
		break;
	    }
	}
        current++;
    }
    func->add( from);
}

// Check whether the declaration starts with a template keyword, and
// if so, use the INLINE_TAG to add the regular expression for an optional
// inline keyword right after the template declaration.
void checkTemplateKeyword( Buffer *decl){
    if( 0 == strncmp( "template", decl->begin(), 8)) {
        const char* p = decl->begin() + 8;
	int nesting = 0;
	while( *p != '\0') {
	    switch( *p) {
	    case '<':
	    case '(':
	        nesting ++;
		break;
	    case ')':
	        nesting --;
		break;
	    case '>':
	        nesting --;
		if ( nesting == 0) {
		    decl->set( p - decl->begin(), INLINE_TAG);
		    return;
		}
		break;
	    }
	    p++;
	}
    }
}


/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleMainComment( const Buffer_list& ) {
    return;
}

void handleComment( const Buffer_list& ) {
    return;
}

void handleClass( const char* classname) {
    if ( global_classname)
        free( global_classname);
    global_classname = strdup( classname);
    errmessage->flush();
    pattern->flush();

    const char* s = classname;
    while ( *s != 0 && *s != '<') 
	s++;
    if ( *s == 0) {
	errmessage->add( "class ");
	errmessage->add( classname);
	errmessage->add( " ");
	pattern->add( "class{ws}+\"");
	pattern->add( classname);
	pattern->add( "\"{classend}");
    } else {
	global_template_params = global_classname;
	while( *global_template_params && *global_template_params != '<')
	    ++global_template_params;
	errmessage->add( "template < ...");
	pattern->add(    "template{ws}*\"<\"{optionalclass}\"");

        int nesting = 0;
	s++;
	while ( nesting >= 0 && *s != 0) {
	    switch ( *s) {
	    case '<':
	        nesting ++;
		errmessage->add( *s);
		pattern->add(    patternText( *s, s[1]));
		break;
	    case '>':
	        nesting --;
		if ( nesting >= 0) {
		    errmessage->add( *s);
		    pattern->add(    patternText( *s, s[1]));
		}
		break;
	    case ',':
	        if ( nesting == 0) {
		    errmessage->add( *s);
		    errmessage->add( " ...");
		    pattern->add(    "\"{ws}*\",\"{optionalclass}\"");
		} else {
		    errmessage->add( *s);
		    pattern->add(    patternText( *s, s[1]));
		}
		break;
	    default:
	        errmessage->add( *s);
		pattern->add(    patternText( *s, s[1]));
		break;
	    }
	    s++;
	}
	if ( nesting >= 0)
	    printErrorMessage( MalformedTemplateParamError);

	errmessage->add( "> class ");
	pattern->add(    "\"{ws}*\">\"{ws}*\"class\"{ws}+\"");

	s = classname;
	while ( *s != 0 && *s != '<' && *s != ' ') {
	    errmessage->add( *s);
	    pattern->add(    *s);
	    s++;
	}
	pattern->add(    "\"{classend}");
    }
    printPattern( *p_out_pattern, *p_out_errmessage, errmessage, pattern);
}

void handleClassEnd( void) {
    global_template_params = 0;
    if ( global_classname)
        free( global_classname);
    global_classname = NULL;
    return;
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
    errmessage->flush();
    int l = strlen( decl);
    while ( l > 0 && decl[ l-1] == ' ') l--;
    if ( decl[ l-1] != ';' )
        printErrorMessage( SemicolonMissingError);
    else 
        l--;
    while ( l > 0 && decl[ l-1] == ' ') l--;
    errmessage->add( decl, l);
    printPattern( *p_out_pattern, *p_out_errmessage, errmessage);
}

void handleNestedType( const char* decl) {
    errmessage->flush();
    int l = strlen( decl);
    while ( l > 0 && decl[ l-1] == ' ') l--;
    errmessage->add( decl, l);
    pattern->flush();
    pattern->add( "typedef{ws}*[^;]+");
    pattern->add( decl, l);
    pattern->add( "{ws}*;");
    // Sorry, not working in full general.
    // pattern->add( "{ws}*;)|(((struct)|(class)){ws}*");
    // pattern->add( decl, l);
    // pattern->add( "{ws}*[:\\{;])");
    printPattern( *p_out_pattern, *p_out_errmessage, errmessage, pattern);
}

void handleMethodDeclaration( const char* decl) {
    translateFunctionArgumentLists( decl, function_decl);
    checkTemplateKeyword( function_decl);
    decl = function_decl->begin();
    errmessage->flush();
    int l = strlen( decl);
    while ( l > 0 && decl[ l-1] == ' ') l--;
    if ( decl[ l-1] != ';' )
        printErrorMessage( SemicolonMissingError);
    else 
        l--;
    while ( l > 0 && decl[ l-1] == ' ') l--;
    errmessage->add( decl, l);
    printFunctionPattern( *p_out_pattern, *p_out_errmessage, errmessage);
}

void handleFunctionDeclaration( const char* decl) {
    translateFunctionArgumentLists( decl, function_decl);
    checkTemplateKeyword( function_decl);
    decl = function_decl->begin();
    errmessage->flush();
    int l = strlen( decl);
    while ( l > 0 && decl[ l-1] == ' ') l--;
    if ( decl[ l-1] != ';' )
        printErrorMessage( SemicolonMissingError);
    else 
        l--;
    while ( l > 0 && decl[ l-1] == ' ') l--;
    errmessage->add( decl, l);
    printFunctionPattern( *p_out_pattern, *p_out_errmessage, errmessage);
}

void handleFunctionTemplateDeclaration( const char* templ, const char* decl) {
    translateFunctionArgumentLists( decl, function_decl);
    decl = function_decl->begin();
    errmessage->flush();
    pattern->flush();
    pattern2->flush();
    errmessage->add( "template < ...");
    pattern->add(    "template{ws}*\"<\"{optionalclass}\"");

    const char* s = templ;
    int nesting = 0;
    while ( *s != 0) {
	    switch ( *s) {
	    case '<':
	        nesting ++;
		errmessage->add( *s);
		pattern->add(    patternText( *s, s[1]));
		break;
	    case '>':
	        nesting --;
		if ( nesting >= 0)
		    errmessage->add( *s);
		    pattern->add(    patternText( *s, s[1]));
		break;
	    case ',':
	        if ( nesting == 0) {
		    errmessage->add( *s);
		    errmessage->add( " ...");
		    pattern->add(    "\"{ws}*\",\"{optionalclass}\"");
		} else {
		    errmessage->add( *s);
		    pattern->add(    patternText( *s, s[1]));
		}
		break;
	    default:
	        errmessage->add( *s);
		pattern->add(    patternText( *s, s[1]));
		break;
	    }
	    s++;
    }
    if ( nesting != 0)
        printErrorMessage( MalformedTemplateParamError);
    errmessage->add( "> ... ");
    pattern->add(    "\"{ws}*\">\"{ws}*{optionalkeywords}{ws}*\"");

    int l = strlen( decl);
    while ( l > 0 && decl[ l-1] == ' ') l--;
    if ( decl[ l-1] != ';' )
        printErrorMessage( SemicolonMissingError);
    else 
        l--;
    while ( l > 0 && decl[ l-1] == ' ') l--;
    errmessage->add( decl, l);
    nesting = 0;
    bool default_params = false;
    for ( int i=0; i<l; i++) {
        switch ( decl[ i]) {
	case '(':
	    nesting ++;
	    break;
	case ')':
	    nesting --;
	    if ( nesting == 0) {
	        default_params = true;
		pattern->add( "\"{ws}*\"");
		pattern2->add( * pattern);
		pattern2->add( ",\".*\"");
	    }
	    break;
	}
	pattern->add( patternText( decl[i], decl[i+1]));
	if ( default_params)
	    pattern2->add( patternText( decl[i], decl[i+1]));
    }
    pattern->add( '"');
    if ( default_params)
        pattern2->add( '"');
    printFunctionPattern( *p_out_pattern, 
			  *p_out_errmessage, 
			  errmessage, 
			  pattern,
			  default_params ? pattern2 : 0);
}


/* main */
/* ==== */
 
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
        cerr << "Usage: cc_build_checker [<options>] "
	        "<scanner-outfile> <messages-header-file> [<infile>]" << endl;
        cerr << "       -trace       sets the `yydebug' variable of bison"
             << endl;
        exit( 1);
    }

    p_out_pattern = new ofstream( parameters[0]);
    if ( ! (*p_out_pattern)) {
        cerr << "error: cannot open file `" << parameters[0] << "' to write."
	     << endl;
	exit( 1);
    }
    p_out_errmessage = new ofstream( parameters[1]);
    if ( ! (*p_out_errmessage)) {
        cerr << "error: cannot open file `" << parameters[1] << "' to write."
	     << endl;
	exit( 1);
    }
    error_filename = parameters[ 1];
    if ( nParameters > 2)
        file_name = parameters[2];
    printProlog( *p_out_pattern, *p_out_errmessage);

    if ( nParameters > 2) {
        for ( i = 2; i < nParameters; i++) {
	    FILE* in;
            if ( (in = fopen( parameters[i], "r")) == NULL) {
	        fprintf( stderr, 
			 "error: cannot open file `%s' for reading.\n",
                         parameters[i]);
                exit( 1);
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

    printEpilog( *p_out_pattern, *p_out_errmessage);
    if ( ! (*p_out_pattern)) {
        cerr << "error: cannot read file `" << parameters[0] << "'."
	     << endl;
	exit( 1);
    }
    delete p_out_pattern;
    if ( ! (*p_out_errmessage)) {
        cerr << "error: cannot read file `" << parameters[1] << "'."
	     << endl;
	exit( 1);
    }
    delete p_out_errmessage;

    return 0;
}

// EOF //

