/**************************************************************************
 
  cc_extract_html.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
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
#include <fstream.h>
#include <strstream.h>
#include <ctype.h>
#include <assert.h>

#include <map.h>
// Repair what the STL bool definitions made wrong.
#ifdef bool
#undef bool
#endif
#ifdef true
#undef true
#endif
#ifdef false
#undef false
#endif

#include <database.h>
#include <html_config.h>

/* Constants:                        */
/* ================================= */
const char* prog_name    = "cc_extract_html";
const char* prog_release = "$Revision$";

/* Configurable command line options */
/* ================================= */
typedef char Switch;
 
#define NO_SWITCH    0
#define MINUS_SWITCH 1
#define PLUS_SWITCH  2
 
Switch  trace_switch  = NO_SWITCH;
Switch  line_switch   = NO_SWITCH;

Switch  config_switch = NO_SWITCH;
Switch  warn_switch   = NO_SWITCH;
Switch  macro_switch  = NO_SWITCH;


/* Two queues to manage footnotes. */
/* =============================== */
int  main_footnote_counter = 0;
int  class_footnote_counter = 0;
int* footnote_counter = NULL;

PQueue<char*>  main_footnotes;
PQueue<char*>  class_footnotes;
PQueue<char*>* footnotes = NULL;

void insertFootnote( char* s) {
    if (footnotes == NULL)
        return;
    footnotes->append(s);
}

// increments counter and returns current value.
int nextFootnoteCounter() {
    if (footnotes == NULL)
        return 0;
    (*footnote_counter)++;
    return *footnote_counter;
}

// format footnote reference and hyperlink based on actual counter 
// into a Buffer.
Buffer* formattedFootnoteNumber() {
    Buffer* buf = new Buffer;
    if (footnotes == NULL)
        buf->add( "[???]");
    else {
        char* tmp = new char[10];
        ostrstream num( tmp, 10);
        num << *footnote_counter << '\0';
        buf->add( "[<A HREF=\"#Footnote_");
        buf->add( tmp);
        buf->add( "\">^");
        buf->add( tmp);
        buf->add( "</A>]");
        delete[] tmp;
    }
    return buf;
}

// prints footnote reference and hyperlink based on actual counter.
void printFootnoteCounter( ostream& out) {
    if (footnotes == NULL)
        out << "[???]";
    else
        out << "[<A HREF=\"#Footnote_" << *footnote_counter << "\">^"
	    << *footnote_counter << "</A>]";
}

// prints footnotes and resets counter.
void printFootnotes( ostream& out) {
    if (footnotes == NULL)
        return;
    if (footnotes->isEmpty())
        return;
    int i = 1;
    out << endl << endl << "<P><HR><H3>Footnotes</H3><P>" << endl << endl;
    out << "<OL>" << endl;
    while( ! footnotes->isEmpty()) {
        out << "<LI><A NAME=\"Footnote_" << i << "\">" 
	    << footnotes->front() << endl;
        delete[] footnotes->front();
        footnotes->remove();
        i++;
    }
    out << "</OL>" << endl;
    *footnote_counter = 0;
}

/* A map for user defined macros     */
/* ================================= */
struct string_less {
    bool operator()(const char* x, const char* y) const { 
        return strcmp(x,y) < 0; 
    }
};

map< const char*, const char*, string_less> user_macros;

void insertMacro( const char* key, const char* value) {
    if ( macro_switch) {
        cerr << "Macro definition: macro `" << key << "' = `" 
	     << value << "'." << endl;
    }
    user_macros[ key] = value;
}

const char* fetchMacro( const char* key) {
    map< const char*, const char*, string_less>::iterator 
        i = user_macros.find(key);
    if ( macro_switch) {
        cerr << "Macro query: macro `" << key << "' = `";
        if (i == user_macros.end())
	    cerr << "<undefined>'." << endl;
        else
	    cerr << (*i).second << "'." << endl;
    }
    if (i == user_macros.end())
        return 0;
    return (*i).second;
}

// The string is already converted to HTML, i.e. $...$ --> <MATH>...</MATH>.
char* extractRCS( char* s) {
    char* old = s;
    while( *s != '>') {
        ADT_Assert( *s);
        ++s;
    }
    ++s;
    while( *s == ' ') {
        ADT_Assert( *s);
        ++s;
    }
    char* p = s;
    while( *s != '<') {
        ADT_Assert( *s);
        ++s;
    }
    --s;
    while( *s == ' ') {
        ADT_Assert( *s);
        --s;
    }
    ++s;
    *s = '\0';
    char* ret = newstr( p);
    delete[] old;
    return ret;
}

char* extractRCSDate( char* s) {
    char* p = extractRCS( s);
    char *q = p + strlen( p) - 1;
    while( *q != ' ') {
        ADT_Assert( *q);
        --q;
    }
    while( *q == ' ') {
        ADT_Assert( *q);
        --q;
    }
    ++q;
    *q = '\0';
    char* ret = newstr( p);
    delete[] p;
    return ret;
}

/* Configurable command line options */
/* ================================= */

// This factor is multiplied to the actual width of an C++ declaration
// right before the test for multiple lines formatting occurs.
// A greater value forces declarations to be printed in multiple lines.
double stretch_factor = 1.6;  // The command line option -strech multiplies
                              // a factor to this builtin default.

// Path for the HTML conversion tools for the default configuration files.
// This path will be compiled into the cc_extract_html program. It is set 
// in the Makefile. The same variable has to be configured in the 
// cc_manual_to_html script.
#ifndef HTML_DEFAULT_PATH
#define HTML_DEFAULT_PATH   NULL
#endif
const char* config_path    = HTML_DEFAULT_PATH;

// Manual date, release, title, and author used in filtering the config files.
const char* manual_date    = NULL;
const char* manual_release = NULL;
const char* manual_title   = NULL;
const char* manual_author  = NULL;

// Directory for the temporary files. A default is given.
// This directory is the output directory for this program. Usually, 
// the script cc_manual_to_html makes a own tmp directory and removes 
// it afterwards. The path must terminate with a slash.
const char* tmp_path       = "/usr/tmp/";


/* Customization tags for the style */
/* ================================ */
bool tag_chapter_author;
bool tag_replace_prefix;     // not supported yet
bool tag_replace_include;    // not supported yet
bool tag_long_param_layout;  // not supported yet

bool tag_rm_const_ref_pair;
bool tag_rm_eigen_class_name;
bool tag_operator_layout;
bool tag_rm_trailing_const;

bool tag_rm_template;
bool tag_template_inline;    // not supported yet

void tag_defaults() {
    tag_chapter_author         = true;
    tag_replace_prefix         = false;
    tag_replace_include        = false;
    tag_long_param_layout      = false;

    tag_rm_const_ref_pair      = true;
    tag_rm_eigen_class_name    = true;
    tag_operator_layout        = true;
    tag_rm_trailing_const      = true;

    tag_rm_template            = false;
    tag_template_inline        = false;
}

void tag_full_declarations() {
    tag_rm_const_ref_pair      = false;
    tag_rm_eigen_class_name    = false;
    tag_operator_layout        = false;
    tag_rm_trailing_const      = false;
    tag_rm_template            = false;
}

void handleNewCommand( char* idfier, char* body){
    // check the body for true or false
    char* p = body;
    while ( *p && isspace( *p))
        p++;
    bool tag = ( 0 == strncmp( p, "\\ccTrue", 7));


    p = idfier;
    while ( *p && isspace( *p))
        p++;
    char *q = p;
    ++q;
    while ( *q && isalpha( *q))
        q++;
    *q = '\0';
    p = newstr( p);
    delete[] idfier;
    insertMacro( p, body);

    // check the different idfiers for tag definitions.
    if ( 0 == strcmp( p, "\\ccTagChapterAuthor")) 
        tag_chapter_author = tag;
    if ( 0 == strcmp( p, "\\ccTagReplacePrefix")) 
        tag_replace_prefix = tag;
    if ( 0 == strcmp( p, "\\ccTagReplaceInclude")) 
        tag_replace_include = tag;
    if ( 0 == strcmp( p, "\\ccAlternateThreeColumn"))  // compatibility
        tag_long_param_layout = !tag;
    if ( 0 == strcmp( p, "\\ccLongParamLayout")) 
        tag_long_param_layout = tag;
    if ( 0 == strcmp( p, "\\ccTagRmTrailingConst")) 
        tag_rm_trailing_const = tag;
    if ( 0 == strcmp( p, "\\ccTagRmEigenClassName")) 
        tag_rm_eigen_class_name = tag;
    if ( 0 == strcmp( p, "\\ccTagRmConstRefPair")) 
        tag_rm_const_ref_pair = tag;
    if ( 0 == strcmp( p, "\\ccTagOperatorLayout")) 
        tag_operator_layout = tag;
    if ( 0 == strcmp( p, "\\ccTagRmTemplate")) 
        tag_rm_template = tag;
    if ( 0 == strcmp( p, "\\ccTagTemplateInline")) 
        tag_template_inline = tag;
}

void handleTexDef( const char* idfier, char* body){
    // idfier starts with \gdef or \def
    idfier++;  // skip first '\'
    while (*idfier && *idfier != '\\')
        idfier ++;
    handleNewCommand( newstr(idfier), body);
}

/* Sort keys for the index            */
/* ================================== */
// The index is organized in a set of fixed topics. 

// scramble is used to write index entries with scope in front.
// The scope is a class but it should not get hyperlinked for its own.
// Has to be same as in the cc_anchor_header.
const char* index_scramble = "<!scramble>";

// sort_keys are used to sort the index according to different sections.
// A sort_key is followed by a 0 to indicate the section title.
// A sort_key is followed by a 1 to indicate a normal entry.

const char* sort_key_class            = "<!sort B";
const char* sort_key_nested_type      = "<!sort C";
const char* sort_key_struct           = "<!sort D";
const char* sort_key_enum             = "<!sort E";
const char* sort_key_enum_tags        = "<!sort F";
const char* sort_key_typedef          = "<!sort G";
const char* sort_key_variable         = "<!sort H";
const char* sort_key_function         = "<!sort I";
const char* sort_key_member_function  = "<!sort J";

void write_headers_to_index( ostream& out){
    out << sort_key_class       << "0 !><P><LI><B>Classes</B><P>" << endl;
    out << sort_key_nested_type << "0 !><P><LI><B>Nested Types</B><P>" << endl;
    out << sort_key_struct      << "0 !><P><LI><B>Structs</B><P>" << endl;
    out << sort_key_enum        << "0 !><P><LI><B>Enums</B><P>" << endl;
    out << sort_key_enum_tags   << "0 !><P><LI><B>Enum Tags</B><P>" << endl;
    out << sort_key_typedef     << "0 !><P><LI><B>Typedefs</B><P>" << endl;
    out << sort_key_variable    << "0 !><P><LI><B>Global Variables and "
                                   "Consts</B><P>" << endl;
    out << sort_key_function    << "0 !><P><LI><B>Functions</B><P>" << endl;
    //    out << sort_key_member_function  << "0 !><P><LI><B>Member Functions</B><P>"
    //	                             << endl;
}

/* table size and font size constants */
/* ================================== */
const int table_width      = 550;
const int table_first_col  = 25;  // in percent
const int table_second_col = 25;  // in percent
const int table_third_col  = 50;  // in percent
const int table_2c_first_col  = 30;  // in percent
const int table_2c_second_col = 70;  // in percent

const double width_per_character  = 5.5;

/* An empty List as empty comment for global declarations */
/* ====================================================== */
Text empty_comment;

/* for the bibliography */
/* ==================== */
bool first_bibitem = true;

/* Declarations from syntax.y */
/* ========================== */
extern int yydebug;
void yyparse();

/* Declarations from lex.yy   */
/* ========================== */
void init_scanner( FILE* in);
extern char* creationvariable;
extern char* formatted_creationvariable;


/* File and filename handling */
/* ========================== */
const char* in_filename = 0;
char* main_filename = "<cout>";
char* current_filename = "<cout>";
char* chapter_title  = 0;
char* class_filename = 0;
char* class_name = 0;
char* template_class_name = 0;
char* formatted_class_name = 0;
char* formatted_template_class_name = 0;
char* scrambled_template_class_name = 0;

const char* html_suffix       = ".html";
const char* chapter_prefix    = "Chapter_";
const char* anchor_filename   = "cc_anchor_rules";
const char* contents_filename = "contents.html";
const char* bib_filename      = "biblio.html";
const char* index_filename    = "cc_index_body";

const char* cc_anchor_header  = "cc_anchor_header";
const char* cc_anchor_footer  = "cc_anchor_footer";
const char* cc_manual_header  = "cc_manual_header";
const char* cc_manual_footer  = "cc_manual_footer";
const char* cc_toc_header     = "cc_toc_header";
const char* cc_toc_footer     = "cc_toc_footer";
const char* cc_biblio_header  = "cc_biblio_header";
const char* cc_biblio_footer  = "cc_biblio_footer";
const char* cc_index_header   = "cc_index_header";
const char* cc_index_footer   = "cc_index_footer";

ostream* main_stream     = 0;
ostream* class_stream    = 0;
ostream* current_stream  = 0;
ostream* anchor_stream   = 0;
ostream* contents_stream = 0;
ostream* index_stream    = 0;

char* addSuffix( const char* name, const char* suffix) {
    char *fullname = new char[ strlen( name) + strlen( suffix) + 1];
    strcpy( fullname, name);
    strcat( fullname, suffix);
    return fullname;
}

// return the position of a suffix dot `.'. Return the terminating
// zero `\0' position if no dot `.' is in the whole filename.
int  suffixPosition( const char* name) {
    int l = strlen( name);
    int i = l;
    while( i > 0) {
        --i;
	if ( name[i] == '.')
	    return i;
    }
    return l;
}

char* replaceSuffix( const char* name, const char* suffix) {
    int l = suffixPosition( name);
    char *fullname = new char[ l + strlen( suffix) + 1];
    strncpy( fullname, name, l);
    fullname[ l] = '\0';
    strcat( fullname, suffix);
    return fullname;
}

char* addPrefix( const char* prefix, const char* name) {
    char *fullname = new char[ strlen( name) + strlen( prefix) + 1];
    strcpy( fullname, prefix);
    strcat( fullname, name);
    return fullname;
}

/* Auxiliary functions for stream handling */
/* ======================================= */

bool exist_file( const char* name){
    istream* in = new ifstream( name);
    if ( ! *in)
        return false;
    return true;
}

void assert_file_write( ostream& out, const char* name){
    if ( ! out) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot write to file `" 
	     << name << "'." << endl;
	exit(1);
    }
}

bool exist_file_with_path( const char* path, const char* name){
    char* full_name = addPrefix( path, name);
    istream* in = new ifstream( full_name);
    delete[] full_name;
    if ( ! *in)
        return false;
    return true;
}

istream* open_file_for_read( const char* name){
    istream* in = new ifstream( name);
    if ( ! *in) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot open file `" << name
	     << "' for reading." << endl;
	exit(1);
    }
    return in;
}

istream* open_file_with_path_for_read( const char* path, const char* name){
    char* full_name = addPrefix( path, name);
    istream* in = new ifstream( full_name);
    if ( ! *in) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot open file `" << full_name
	     << "' for reading." << endl;
	exit(1);
    }
    delete[] full_name;
    return in;
}

ostream* open_file_with_path_for_write( const char* path, const char* name){
    char* full_name = addPrefix( path, name);
    ostream* out = new ofstream( full_name);
    if ( ! *out) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot open file `" << full_name
	     << "' for writing." << endl;
	exit(1);
    }
    delete[] full_name;
    return out;
}


/* read a file into a buffer   */
/* The name has a trailing '}' */
/* HTML conversion is done simultaneously */
/* =========================== */
const char* html_multi_character( char c);  // prototype from below

Buffer* readFileInBuffer(const char* name) {
    int len = strlen( name);
    const char* p = name + len - 1;
    while( p != name && ( *p == ' ' || *p == '}'))
        p--;
    char* filename = new char[p-name+2];
    strncpy( filename, name, p-name+1);
    filename[p-name+1] = '\0';
    istream* in = open_file_for_read(filename);
    delete[] filename;
    char c;
    Buffer* b = new Buffer;
    while( in->get(c))
        b->add( html_multi_character(c));
    return b;
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
		    is_empty = class_name == NULL;
		    break;
		case 'u':
		    is_empty = main_stream == &cout;
		    break;
		case 'd':
		    is_empty = manual_date == NULL;
		    if ( !is_empty)
		        is_empty = (strlen(manual_date) == 0);
		    break;
		case 'r':
		    is_empty = manual_release == NULL;
		    if ( !is_empty)
		        is_empty = (strlen(manual_release) == 0);
		    break;
		case 't':
		    is_empty = manual_title == NULL;
		    if ( !is_empty)
		        is_empty = (strlen(manual_title) == 0);
		    break;
		case 'a':
		    is_empty = manual_author == NULL;
		    if ( !is_empty)
		        is_empty = (strlen(manual_author) == 0);
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
		    if ( class_name)
		        out << template_class_name;
		    break;
		case 'u':
		    if ( main_stream != &cout)
		        out <<  main_filename;
		    break;
		case 'd':
		    if ( manual_date ) 
		        out << manual_date;
		    break;
		case 'r':
		    if ( manual_release ) 
		        out << manual_release;
		    break;
		case 't':
		    if ( manual_title ) 
		        out << manual_title;
		    break;
		case 'a':
		    if ( manual_author ) 
		        out << manual_author;
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

istream* open_config_file( const char* name){
    if ( config_switch == NO_SWITCH) {
        // check if the file exists in the current directory
        if ( exist_file( name))
	    return (open_file_for_read( name));
    }
    return( open_file_with_path_for_read( config_path, name));
}

void copy_and_filter_config_file( const char* name, ostream& out){
    istream* in  = open_config_file( name);
    filter_config_file( *in, out);
    delete in;
}

void copy_config_file( const char* name){
    istream* in  = open_config_file( name);
    ostream* out = open_file_with_path_for_write( tmp_path, name);
    filter_config_file( *in, *out);
    delete in;
    assert_file_write( *out, name);
    delete out;
}


/* Filter decl. before writing to the index comment */
/* ================================================ */

void filter_for_index_comment( ostream& out, const char* text) {
    while( *text) {
        switch (*text) {
	case '<':
	    out << '(';	  
	    break;
	case '>':
	    out << ')';	  
	    break;
	case '&':
	    out << '_';	  
	    break;
	default:
	    out << *text;
	}
        text++;
    }
}

// Filter characters that might not be allowed in hyperlink anchors.
void filter_for_index_anchor( ostream& out, const char* text) {
    while( *text) {
        switch ( *text) {
	case '<':
	case '(':
	case '[':
	case '{':
	    out << '6';
	    break;
	case '>':
	case ')':
	case ']':
	case '}':
	    out << '9';
	    break;
	case '"':
	case '&':
	case ' ':
	    out << '_';
	    break;
	case ',':
	case '.':
	    out << '+';
	    break;
	default:
	    out << *text;
	}
        text++;
    }
}



/* HTML generating functions */
/* ========================= */

void open_html( ostream& out) {
    istream* in = open_config_file( cc_manual_header);
    filter_config_file( *in, out);
    delete in;
}

void close_html( ostream& out) {
    istream* in = open_config_file( cc_manual_footer);
    filter_config_file( *in, out);
    delete in;
}

bool is_html_multi_character( char c) {
    return c == '"' || c == '&' || c == '<' || c == '>';
}

const char* html_multi_character( char c) {
    static char *s = " ";
    switch ( c) {
    case '"': return "&quot;";
    case '&': return "&amp;";
    case '<': return "&lt;";
    case '>': return "&gt;";
    default:  *s = c;
    }
    return s;
}

void print_ascii_to_html( ostream& out, const char* txt) {
    while( *txt) {
        if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    out << *txt;
	++txt;
    }
}

void print_ascii_len_to_html( ostream& out, const char* txt, int n) {
    while( n) {
        if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    out << *txt;
	++txt;
	n--;
    }
}

// This version eliminates multiple spaces.
void print_ascii_to_html_spc( ostream& out, const char* txt) {
    while( *txt) {
        if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    if ( *txt > ' ' || (*txt > '\0' && txt[1] > ' '))
	        out << *txt;
	++txt;
    }
}

int strlen_ascii_to_html( const char* txt) {
    int len = 0;
    while( *txt) {
        if ( is_html_multi_character( *txt))
	    len += strlen( html_multi_character( *txt));
	else
	    len++;
	++txt;
    }
    return len;
}

char* convert_ascii_to_html( const char* txt) {
    if ( txt == NULL) {
        char *q = new char[1];
	q[0] = '\0';
	return q;
    }
    char* s = new char[ strlen_ascii_to_html( txt) + 1];
    char* p = s;
    while( *txt) {
        if ( is_html_multi_character( *txt)) {
	    const char* q = html_multi_character( *txt);
	    while ( *q)
	        *p++ = *q++;
	} else
	    *p++ = *txt;
	++txt;
    }
    *p = '\0';
    return s;
}

// A scope operator is also appended
char* convert_ascii_to_scrambled_html( const char* txt) {
    char* s = new char[  strlen_ascii_to_html( txt) + 3 
                       + strlen( index_scramble)];
    char* p = s;
    bool tag = true;
    while( *txt) {
        if ( is_html_multi_character( *txt)) {
	    const char* q = html_multi_character( *txt);
	    while ( *q)
	        *p++ = *q++;
	} else
	    *p++ = *txt;
	++txt;
	if ( tag) {
	    const char* q = index_scramble;
	    while ( *q)
	        *p++ = *q++;
	    tag = false;
	}
    }
    *p++ = ':';
    *p++ = ':';
    *p = '\0';
    return s;
}

double estimate_html_size( const char* s) {
    int n = 0;
    while ( *s) {
        if ( *s > ' ' || (*s > '\0' && s[1] > ' '))
	    n++;
        ++s;
    }
    return width_per_character * n; 
}


void three_cols_html_begin( ostream& out, bool big_col1) {
    out << indent << indent << indNewline
	<< "<!3><TABLE BORDER=0 CELLSPACING=2 CELLPADDING=0 WIDTH="
	<< table_width << ">" << indNewline 
	<< "<TR><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_first_col << "% NOWRAP" << ( big_col1 ? " COLSPAN=3>" : ">")
	<< indNewline << "<VAR>" << outdent << indNewline;
}

void three_cols_html_premature_end( ostream& out) {
    out << indent << indNewline	<<"</VAR>" << indNewline << "</TD>";
    out << "</TR>" << indNewline 
	<< "</TABLE><!3>" << outdent << outdent << indNewline;
}

void three_cols_html_second( ostream& out, bool big_col1, bool big_col2) {
    out << indent << indNewline	<<"</VAR>" << indNewline << "</TD>";
    if ( big_col1)
        out << "</TR><TR><TD WIDTH=" << table_first_col << "%></TD>";
    out << "<TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_second_col << "% NOWRAP" << ( big_col2 ? " COLSPAN=2>" : ">")
	<< indNewline << "<VAR>" << outdent << indNewline;
}

void three_cols_html_third( ostream& out, bool big_col2, bool empty_col3) {
    out << indent << indNewline << "</VAR>" << indNewline << "</TD>";
    if ( big_col2 && ! empty_col3)
        out << "</TR><TR><TD WIDTH=" << table_first_col 
	    << "%></TD><TD WIDTH=" << table_second_col << "%></TD>";
    out << "<TD ALIGN=LEFT VALIGN=TOP WIDTH=" << table_third_col << "%>" ;
    out << outdent << indNewline;
}

void three_cols_html_end( ostream& out, bool big_col2, bool empty_col3) {
    out << indent << indNewline;
    if ( ! big_col2 || ! empty_col3)
        out << "</TD>";
    out << "</TR>" << indNewline 
	<< "</TABLE><!3>" << outdent << outdent << indNewline;
}

void two_cols_html_begin( ostream& out) {
    out << indent << indent << indNewline
	<< "<!2><TABLE BORDER=0 CELLSPACING=2 CELLPADDING=0 WIDTH="
	<< table_width << ">" << indNewline 
	<< "<TR><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_2c_first_col << "% NOWRAP COLSPAN=2>"
	<< indNewline << "<VAR>" << outdent << indNewline;
}

void two_cols_html_premature_end( ostream& out) {
    out << indent << indNewline	<<"</VAR>" << indNewline << "</TD>";
    out << "</TR>" << indNewline 
	<< "</TABLE><!2>" << outdent << outdent << indNewline;
}

void two_cols_html_second( ostream& out, bool empty_col2) {
    out << indent << indNewline << "</VAR>" << indNewline << "</TD></TR>";
    if ( ! empty_col2)
        out << "<TR><TD WIDTH=" << table_2c_first_col 
	    << "%></TD><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	    << table_2c_second_col << "%>";
    out << outdent << indNewline;
}

void two_cols_html_end( ostream& out, bool empty_col2) {
    out << indent << indNewline;
    if ( ! empty_col2)    
        out << "</TD></TR>" << indNewline;
    out	<< "</TABLE><!2>" << outdent << outdent << indNewline;
}

int text_block_length( const Text& T) {
    int l = 0;
    InListFIter< TextToken> words( (Text&)T);
    ForAll( words) {
        if ( words->isSpace)
	    ++l;
	else
	    l += (*words).len;
    }
    return l;
}

char* text_block_to_string( const Text& T) {
    char* string = new char[ text_block_length( T) + 1];
    string[0] = '\0';
    InListFIter< TextToken> words( (Text&)T);
    ForAll( words) {
        if ( words->isSpace)
	    strcat( string, " ");
	else
	    strcat( string, (*words).string);
    }
    return string;
}

bool is_text_block_empty( const Text& T) {
    InListFIter< TextToken> words( (Text&)T);
    ForAll( words) {
        if ( !words->isSpace)
	    return false;
    }
    return true;
}

int print_html_text_block( ostream &out, const Text& T,
			   bool leadingParagraph = false) {
    InListFIter< TextToken> words( (Text&)T);
    int  width = MaxTextWidth - indentation_number() - 3;
    int  w     = width;
    int  state = 0;     // 0 = start, 1 = after token, 2 = spaces
                        // 3 = one newline (and spaces), 4 = newlines
    ForAll( words) {
	switch ( state) {
	case 0:
	    if ( !words->isSpace) {
	        if ( leadingParagraph) {
		    out << indNewline << "<P>\n" << indNewline;
		}
		out << words->string;
		w -= words->len;
		state = 1;
	    }
	    break;
	case 1:
	    if ( !words->isSpace || words->len > 0) {
		if ( !words->isSpace) {
		    if ((words->len > w) && (w != width)) {
			w = width;
			out << indNewline;
		    }
		    out << words->string;
		    w -= words->len;
		} else {
		    if ( words->string[0] == '\n')
			state = 3;
		    else
			state = 2;
		}
	    }
	    break;
	case 2:
	    if ( !words->isSpace || words->len > 0) {
		if ( !words->isSpace) {
		    if ((words->len >= w) && (w != width)) {
			w = width;
			out << indNewline;
		    } else {
			out << ' ';
			w--;
		    }
		    out << words->string;
		    w -= words->len;
		    state = 1;
		} else {
		    if ( words->string[0] == '\n')
			state = 3;
		}
	    }
	    break;
	case 3:
	    if ( !words->isSpace || words->len > 0) {
		if ( !words->isSpace) {
		    if ((words->len >= w) && (w != width)) {
			w = width;
			out << indNewline;
		    } else {
			out << ' ';
			w--;
		    }
		    out << words->string;
		    w -= words->len;
		    state = 1;
		} else {
		    if ( words->string[0] == '\n')
			state = 4;
		}
	    }
	    break;
	case 4:
	    if ( !words->isSpace || words->len > 0) {
		if ( !words->isSpace) {
		    out << indNewline << "<P>\n" << indNewline;
		    w = width;
		    out << words->string;
		    w -= words->len;
		    state = 1;
		}
	    }
	    break;
	}
    }
    return state;
}


/* template declarations for functions and others */
/* ============================================== */
/* It depends on the layout (2 or 3 columns)      */

const char* handle_template_layout( ostream &out, 
				    const char* decl, 
				    bool three_col) {
    const char* new_pos = decl;
    const char* p = decl;
    while ( *p && isspace( *p))
        p++;
    if ( strncmp( "template", p, 8) == 0) {
        p += 8;
	while ( *p && isspace( *p))
	    p++;
	if ( *p == '<') {
	    p++;
	    // Here we know that it is a template declaration. Now we
            // look for the end of the template params using a nesting counter.
	    int nesting = 1;
	    while ( *p && nesting > 0) {
	        if ( *p == '<') 
		    nesting++;
		else if ( *p == '>')
		    nesting --;
		p++;
	    }
	    new_pos = p;
	    if ( ! tag_rm_template) {
		// Print the template declaration.
		if ( three_col)
		    three_cols_html_begin( out, true);
		else
		    two_cols_html_begin( out);

		print_ascii_len_to_html( out, decl, p - decl);

		if ( three_col)
		    three_cols_html_premature_end( out);
		else
		    two_cols_html_premature_end( out);
	    }
	}
    }
    return new_pos;
}

/* Formatting functions that work like the cc_manual.sty */
/* ===================================================== */

// This function separates the return value from the scope, the function name,
// the parameter list (without paranthesis) and the rest of the signature.
// Therefore, the string in signature will be parsed from the end
// looking for the first ')' and then for the matching '(', using
// a nesting count for '(', ')', '{', '}', '<', '>'. Then, the next token is 
// assumed to be the function name. After that, a scope can be a list of 
// '::' and identifiers. The rest is assumed to be the return value. All parts
// might be empty. The parantheses for the parameters are not allowed to be
// omitted. The function name token might be an operator !!
// For the special case of enumerators the last parameter enum_decl
// changes the parsing to look for braces `{}' instead of parantheses `()'.
void split_function_declaration( const char* signature, 
				 char*& return_value, 
				 char*& scope, 
				 char*& function_name, 
				 char*& parameter_list,
				 const char*& rest,
				 bool  enum_decl = false) {
    char opening = '(';
    char closing = ')';
    if ( enum_decl) {
        opening = '{';
	closing = '}';
    }
    const char* s = signature + strlen( signature) - 1;

    // Note that the string processing pointer s can point one character
    // in front of signature (where it should never be dereferenced).
    // This is checked using the s_end pointer.
    const char* s_end = signature - 1;

    // skip the rest
    while ( s != s_end && *s != closing)
        --s;
    rest = s + 1;
    while ( *rest > 0 && *rest <= ' ')
        ++rest;
    if ( s != s_end)
        --s;
    while ( s != s_end && *s <= ' ')
        --s;

    // scan the parameter list
    const char* q = s;
    int nesting = 0;
    while ( s != s_end && ( nesting || *s != opening)) {
        switch ( *s) {
	case ')':
	case '}':
	case '>':
	    ++nesting;
	    break;
	case '(':
	case '{':
	case '<':
	    --nesting;
	    break;

	}
        --s;
    }
    const char* p = s + 1;
    while ( *p > 0 && *p <= ' ')
        ++p;
    if ( q - p >= 0) {
        parameter_list = new char[ q - p + 2];
        strncpy( parameter_list, p, q - p + 1);
        parameter_list[ q - p + 1] = '\0';
    } else 
        parameter_list = 0;
    
    if ( s != s_end) // read over the '('
        --s;
    while ( s != s_end && *s <= ' ')
        --s;
    if ( s == s_end || nesting) {
        printErrorMessage( MalformedFunctionDeclaration);
	exit( 1);
    }
    q = s;

    // scan function name
    // skip possible operator symbols and white spaces
    while ( s != s_end && ( ! isalnum( *s)) && ( *s != '_'))
        --s;
    // parse first identifier for the function_name
    while ( s != s_end && ( isalnum( *s) || *s == '_'))
        --s;
    // check the possibilty that the function is a cast operator
    // like `operator int', or an idfier as operator like `operator new'
    // so check whether the next identifier is equal to `operator'
    p = s;
    while ( s != s_end && *s <= ' ')
        --s;
    while ( s != s_end && ( isalnum( *s) || *s == '_'))
        --s;
    if ( strncmp( s + 1, "operator", 8) != 0)
        // it's not the cast operator, restore old positiom
        s = p;
    if ( q - s) {
        function_name = new char[ q - s + 1];
        strncpy( function_name, s + 1, q - s);
        function_name[ q - s] = '\0';
    } else 
        function_name = 0;
    
    // check for a scope operator
    while ( s != s_end && *s <= ' ')
        --s;
    q = s;
    while ( (s - s_end) > 0 && *s == ':' && *(s-1) == ':') {
        s -= 2;
	while ( s != s_end && *s <= ' ')
	    --s;
	while ( s != s_end && ( isalnum( *s) || *s == '_'))
	    --s;
    }
    if ( q - s > 0) {
        scope = new char[ q - s + 1];
        strncpy( scope, s + 1, q - s);
        scope[ q - s] = '\0';
    } else 
        scope = 0;

    // The rest is the return type
    while ( s != s_end && *s <= ' ')
        --s;
    q = signature;
    while ( *q && *q <= ' ')
        ++q;
    if ( s - q >= 0) {
        return_value = new char[ s - q + 2];
        strncpy( return_value, q, s - q + 1);
        return_value[ s - q + 1] = '\0';
    } else 
        return_value = 0;
}

// This function separates a variable declaration around the variable name.
// The rest might contain `)', array subscripts or an assignment. 
// The scope might be empty, the return_value should not be empty. 
// No fancy tricks like operators.
void split_variable_declaration( const char* signature, 
				 char*& return_value, 
				 char*& scope, 
				 char*& variable_name, 
				 char*& rest) {
    const char* s = signature + strlen( signature) - 1;

    // Note that the string processing pointer s can point one character
    // in front of signature (where it should never be dereferenced).
    // This is checked using the s_end pointer.
    const char* s_end = signature - 1;

    // skip the `;'
    while ( s != s_end && *s != ';')
        --s;
    if ( s != s_end)
        --s;

    // skip the rest
    const char* q = s;

    // look out for an assignment (including nesting)
    int nesting = 0;
    while ( s != s_end && ( nesting || *s != '=')) {
        if ( *s == '}' || *s == ')' || *s == ']' || *s == '>')
	    ++nesting;
        if ( *s == '{' || *s == '(' || *s == '[' || *s == '<')
	    --nesting;
        --s;
    }
    if ( s != s_end) {
        // assignment found
        --s;
    } else {
        // not found, restore old position
        s = q;
    }

    // skip possible operator symbols and white spaces
    while ( s != s_end && ( ! isalnum( *s)) && ( *s != '_'))
        --s;
    const char* p = s + 1;
    while ( *p > 0 && *p <= ' ')
        ++p;
    if ( q - p >= 0) {
        rest = new char[ q - p + 2];
        strncpy( rest, p, q - p + 1);
        rest[ q - p + 1] = '\0';
    } else
        rest = 0;
    
    if ( s == s_end) {
        printErrorMessage( MalformedFunctionDeclaration);
	exit( 1);
    }
    q = s;

    // scan function name
    // parse first identifier for the function_name
    while ( s != s_end && ( isalnum( *s) || *s == '_'))
        --s;
    if ( q - s) {
        variable_name = new char[ q - s + 1];
        strncpy( variable_name, s + 1, q - s);
        variable_name[ q - s] = '\0';
    } else 
        variable_name = 0;
    
    // check for a scope operator
    while ( s != s_end && *s <= ' ')
        --s;
    q = s;
    while ( (s - s_end) > 0 && *s == ':' && *(s-1) == ':') {
        s -= 2;
	while ( s != s_end && *s <= ' ')
	    --s;
	while ( s != s_end && ( isalnum( *s) || *s == '_'))
	    --s;
    }
    if ( q - s > 0) {
        scope = new char[ q - s + 1];
        strncpy( scope, s + 1, q - s);
        scope[ q - s] = '\0';
    } else 
        scope = 0;

    // The rest is the return type
    while ( s != s_end && *s <= ' ')
        --s;
    q = signature;
    while ( *q && *q <= ' ')
        ++q;
    if ( s - q >= 0) {
        return_value = new char[ s - q + 2];
        strncpy( return_value, q, s - q + 1);
        return_value[ s - q + 1] = '\0';
    } else 
        return_value = 0;
}

void remove_const_ref_pair( char* s) {
    if ( ! tag_rm_const_ref_pair)
        return;
    int state = 0; // 0 = no const, 1 = c, 2 = co, 3 = con, 4 = cons, 5 = const
    int nesting = 0;
    char* q;       // position of recently found const
    int   nest;    // nesting level of recently found const
    while( *s) {
        switch( state) {
	case 0:
	    if ( *s == '(' || *s == '<' || *s == '[' || *s == '{')
	        ++nesting;
	    if ( *s == ')' || *s == '>' || *s == ']' || *s == '}')
	        --nesting;
	    if ( *s == 'c') {
	        state = 1;
		nest = nesting;
		q = s;
	    }
	    ++s;
	    break;
	case 1:
	    if ( *s == 'o') {
	        ++s;
		state = 2;
	    } else
	        state = 0;
	    break;
	case 2:
	    if ( *s == 'n') {
	        ++s;
		state = 3;
	    } else
	        state = 0;
	    break;
	case 3:
	    if ( *s == 's') {
	        ++s;
		state = 4;
	    } else
	        state = 0;
	    break;
	case 4:
	    if ( *s == 't') {
	        ++s;
		state = 5;
	    } else
	        state = 0;
	    break;
	case 5:
	    if ( *s == '(' || *s == '<' || *s == '[' || *s == '{')
	        ++nesting;
	    if ( *s == ')' || *s == '>' || *s == ']' || *s == '}')
	        --nesting;
	    if ( *s == '&' && nest == nesting) {
	        q[0] = ' ';
	        q[1] = ' ';
	        q[2] = ' ';
	        q[3] = ' ';
	        q[4] = ' ';
		s[0] = ' ';
		return;
	    }
	    ++s;
	    break;
	}
    }
    return;
}

void remove_own_classname( char* s, const char* classname) {
    if ( ! tag_rm_eigen_class_name)
        return;
    if ( ! classname)
        return;
    // An easy algorithm (quadratic runtime worst case).
    // The classname has to match exactly (and is removed at most once).
    int l = strlen( s) - strlen( classname);
    int l2 = strlen( classname);
    int i;
    for ( i = 0; i <= l; i++) {
        if ( strncmp( s, classname, l2) == 0)
	    break;
	s++;
    }
    if ( i <= l) {
        // Classname found. Remove it.
	for( ; l2 > 0; l2--)
	    *s++ = ' ';
    }
}

// This function separates each parameter in a parameter list.
// It replaces therefore the separating commatas by `\0'.
// The returnvalue states the number of parameters found.
// The parsing process takes nesting in account.
int separate_parameter_list( char* s) {
    if ( ! s) 
        return 0;
    char* s_start = s - 1;
    int n = 1;
    int nesting = 0;
    while( *s) {
        if ( *s == '(' || *s == '<' || *s == '[' || *s == '{')
	    ++nesting;
	if ( *s == ')' || *s == '>' || *s == ']' || *s == '}')
	    --nesting;
	if ( *s == ',' && nesting == 0) {
	    ++n;
	    // remove trailing spaces
	    *s = ' ';
	    char* q = s - 1;
	    while ( q != s_start && *q && *q <= ' ')
	        --q;
	    ++q;
	    *q = '\0';
	}
	++s;
    }
    if ( nesting ) {
        printErrorMessage( MalformedFunctionDeclaration);
	exit( 1);
    }
    return n;
}

// Compute the operator layout. Assume reasonable values in all parameters.
// Actually, exp_size just counts all parameters and operator characters.
// Only exceptional cases has to recompute exp_size. The boolean 
// ignore_params can be set true in which case the parameters will not
// be printed, even not the creationvariable for methods. (Used
// for the new and the delete operators.)
// The parantheses operator has to format its own parantheses and
// the first commata between the first and second parameter (if
// more than one parameter exists).
// As an exceptional case for ++ and -- the number of parameters 
// is reduced from two to one to hide the second int parameter
// neccessary to indicate that the postincrement is meant. (Hack!! ;--> )
// Returns: false if failed, true for success.
bool format_operators( int n, int& modifiable_n, 
		       char* op,
		       const char*& praefix, 
		       const char*& infix, 
		       const char*& postfix,
		       double& exp_size, 
		       bool& ignore_params) {
    // Assume that the syntax is not malformed (error messages are
    // provided elsewhere, namely compiler or cc_manual.sty).
    // Non matching operators are printed in functional notation.
    if ( ! tag_operator_layout)
        return false;
    if ( *op <= ' ')
        return false;
    if ( op[1] <= ' ') {  // one character operators
        if ( n == 1) {        // with one parameter
	    switch ( *op) {
	    case '~':
	    case '!':
	    case '-':
	    case '+':
	    case '&':
	    case '*':
	        praefix = op;
		break;
	    default:
		return false;
	    }
	} else if ( n == 2) {  // with two parameters
	    switch ( *op) {
	    case '*':
	    case '/':
	    case '%':
	    case '+':
	    case '-':
	    case '<':
	    case '>':
	    case '&':
	    case '^':
	    case '|':
	    case '=':
	    case ',':
	        infix = op;
		break;
	    default:
		return false;
	    }
	} else
	    return false;
    } else if ( op[2] <= ' ') {  // two character operators
        if ( n == 1) {        // with one parameter
	    switch ( *op) {
	    case '-':
	        if ( op[1] == '>')
		    postfix = op;
	        else if ( op[1] == '-')
		    praefix = op;
		else
		    return false;
		break;
	    case '(':
	        if ( op[1] == ')') {
		    infix   = op;
		} else
		    return false;
		break;
	    case '+':
	        if ( op[1] == *op)
		    praefix = op;
		else
		    return false;
		break;
	    default:
		return false;
	    }
	} else if ( n == 2) {  // with two parameters
	    switch ( *op) {
	    case '[':
	        if ( op[1] == ']') {
		    infix   = "[";
		    postfix = "]";
		}
		else
		    return false;
		break;
	    case '(':
	        if ( op[1] == ')') {
		    infix   = "(";
		    postfix = ")";
		}
		else
		    return false;
		break;
	    case '+':
	    case '-':
	        if ( op[1] == *op) {
		    postfix = op;
		    modifiable_n--;
		    exp_size -= estimate_html_size( " int");
		}
	        else if ( op[1] == '=')
		    infix = op;
		else
		    return false;
		break;
	    case '>':
	    case '<':
	    case '&':
	    case '|':
	        if ( op[1] == *op || op[1] == '=')
		    infix = op;
		else
		    return false;
		break;
	    case '=':
	    case '!':
	    case '*':
	    case '/':
	    case '%':
	    case '^':
	        if ( op[1] == '=')
		    infix = op;
		else
		    return false;
		break;
	    default:
		return false;
	    }
	} else
	    return false;
    } else {  // three or more character operators
        if ( strcmp( op, "->*") == 0)
	    postfix = op;
	else if ( strcmp( op, "<<=") == 0)
	    infix   = op;
	else if ( strcmp( op, ">>=") == 0)
	    infix   = op;
	else if ( strcmp( op, "new") == 0) {
	    ignore_params = true;
	    exp_size  = estimate_html_size( " new ");
	    exp_size += estimate_html_size( template_class_name);
	    praefix   = "new";
	    infix     = template_class_name;
	} else if ( strcmp( op, "delete") == 0) {
	    ignore_params = true;
	    exp_size  = estimate_html_size( " delete void*");
	    praefix   = "delete void*";
	} else if ( strcmp( op, "delete[]") == 0) {
	    ignore_params = true;
	    exp_size  = estimate_html_size( " delete[] void*");
	    praefix   = "delete[] void*";
	} else
	    return false;
    }
    return true;
}


void print_rest( ostream& out, const char* txt){
    out << ' ';
    while( *txt && *txt != ';') {
        if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    if ( *txt > ' ' || (*txt > '\0' && txt[1] > ' '))
	        out << *txt;
	++txt;
    }
}

void format_function( bool method, const char* signature, const Text& T) {
    char* return_value;
    char* scope;
    char* function_name; 
    char* parameter_list; 
    const char* rest;

    bool  normal_operator     = false;  // either operator ...
    bool  conversion_operator = false;  // ... or conversion operator 
                                        // or (normal) function
    char* op_symbols;

    while( *signature && *signature == ' ')
        signature++;

    split_function_declaration( signature, 
				return_value,
				scope,
				function_name,
				parameter_list,
				rest);

    // check function_name for operator
    if (     strncmp( function_name, "operator", 8) == 0
	  && ! isalnum(function_name[8]) 
	  && function_name[8] != '_') {
	op_symbols = function_name + 8;
	while ( *op_symbols && *op_symbols <= ' ')
	    ++op_symbols;
        normal_operator     = return_value;
	conversion_operator = ! return_value;
    }

    double exp_size_ret = 0.0;
    double exp_size = 0.0;
    if ( return_value) {
        remove_const_ref_pair( return_value);
	exp_size_ret += estimate_html_size( return_value);
    }
    if ( conversion_operator)
	exp_size_ret += estimate_html_size( op_symbols);
 
    three_cols_html_begin( *current_stream, 
			   exp_size_ret > table_width*table_first_col/100.0);
    // ---------
    // index
    if ( !method) {
	char* formatted_return   = convert_ascii_to_html( return_value);
	char* formatted_scope    = convert_ascii_to_html( scope);
	char* formatted_function = convert_ascii_to_html( function_name);
	char* formatted_params   = convert_ascii_to_html( parameter_list);
	char* formatted_rest     = convert_ascii_to_html( rest);
	// Make only the function name a hyperlink. Types could be substituted
	// according the rules.
	*index_stream << sort_key_function << '1';
	filter_for_index_comment( *index_stream, function_name);
	*index_stream << ' ';
	filter_for_index_comment( *index_stream, signature);
	*index_stream << "!><UL><LI><VAR>" 
		      << formatted_return << ' ' << formatted_scope
		      << "<A HREF=\"" << current_filename << "#Function_";
	filter_for_index_anchor( *index_stream, signature);
	*index_stream << "\">" << formatted_function << "</A>"
		      << "( " << formatted_params << ")" << formatted_rest
		      << "</VAR></UL>" << endl;
	*current_stream << "<A NAME=\"Function_";
	filter_for_index_anchor( *current_stream, signature);
	*current_stream << "\"></A>" << endl;
	delete[] formatted_return;
	delete[] formatted_scope;
	delete[] formatted_function;
	delete[] formatted_params;
	delete[] formatted_rest;
    }
    // end index
    // ----------
   
    if ( return_value)
        print_ascii_to_html_spc( *current_stream, return_value);
    if ( conversion_operator)
        print_ascii_to_html_spc( *current_stream, op_symbols);

    // handle function body or operation signature
    // first, estimate size
    if ( method)
        exp_size = estimate_html_size( creationvariable) + width_per_character;
    int n = 0;
    if ( conversion_operator) {
        exp_size += estimate_html_size( op_symbols);
    } else  if ( parameter_list) {
        n = separate_parameter_list( parameter_list);
	char* p = parameter_list;
	int m = n;
	while ( m--) {
	    remove_const_ref_pair( p);
	    remove_own_classname( p, template_class_name);
	    exp_size += estimate_html_size( p) + width_per_character;
	    p += strlen( p) + 1;  // skip to next parameter
	}
    }
    if (!tag_rm_trailing_const)
        exp_size += estimate_html_size( rest);

    bool failed = false;
    if ( normal_operator) {
	exp_size += estimate_html_size( op_symbols);
	bool ignore_params = false;  // exception for new and delete operators
	const char* praefix = "";
	const char* infix   = "";
	const char* postfix = "";
	failed = ! format_operators( ( method ? n + 1 : n), n,
				     op_symbols,
				     praefix, infix, postfix,
				     exp_size, ignore_params);
	if ( ! failed) {
	    // print the operator
	    three_cols_html_second(
		    *current_stream,
		    exp_size_ret > table_width * table_first_col / 100.0,
		    exp_size > table_width * table_second_col / 100.0
		);
	    print_ascii_to_html_spc( *current_stream, praefix);
	    *current_stream << " ";
	    char* p = parameter_list;
	    if ( ! ignore_params) {
	        if ( method)
		    print_ascii_to_html_spc(*current_stream, creationvariable);
		else if (n) {
		    --n;
		    print_ascii_to_html_spc(*current_stream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		}
	    }
	    *current_stream << " ";
	    print_ascii_to_html_spc( *current_stream, infix);
	    *current_stream << " ";
	    if ( ! ignore_params && n) {
	        while ( n--) {
		    print_ascii_to_html_spc(*current_stream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		    if ( n)
		        *current_stream << ", ";
		}
	    }
	    print_ascii_to_html_spc( *current_stream, postfix);
	}
    }
    if ( ! normal_operator || failed) {
	if ( scope)
	    exp_size += estimate_html_size( scope);
	exp_size += estimate_html_size( function_name);
	exp_size += 3 * width_per_character;  // for parameter list parantheses

	// then, do the printing
	three_cols_html_second(
		*current_stream,
		exp_size_ret > table_width * table_first_col / 100.0,
		exp_size > table_width * table_second_col / 100.0
	    );
	if ( conversion_operator) {
	    print_ascii_to_html_spc( *current_stream, op_symbols);
	    *current_stream << " ( ";
	    print_ascii_to_html_spc( *current_stream, creationvariable);
	    *current_stream << ")";
	} else {
	    double dd_width = table_width * ( 1.0 - table_first_col / 100.0);
	    dd_width /= stretch_factor;
	    if ( exp_size > dd_width && parameter_list)
	        *current_stream<<"<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
		                 "<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
			       << indNewline;
	    if ( method) {
		print_ascii_to_html_spc( *current_stream, creationvariable);
		*current_stream << '.';
	    }
	    if ( scope)
		print_ascii_to_html_spc( *current_stream, scope);
	    print_ascii_to_html_spc( *current_stream, function_name);
	    if ( parameter_list) {
		*current_stream << " ( ";
		if ( exp_size > dd_width)
		    *current_stream <<"</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
				    << indNewline;
		char* p = parameter_list;
		while ( n--) {
		    print_ascii_to_html_spc( *current_stream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		    if ( n) {
			*current_stream << ", ";
			if ( exp_size > dd_width)
			    *current_stream << "<BR>" << indNewline;
		    }
		}
		*current_stream << ")";
		if ( exp_size > dd_width)
		    *current_stream << "</TD></TR></TABLE>" << indNewline;
	    } else
		*current_stream << " ()";
	}
    }
    if (!tag_rm_trailing_const)
        print_rest( *current_stream, rest);
    bool is_empty_comment = is_text_block_empty( T);
    three_cols_html_third( *current_stream, 
			   exp_size > table_width * table_second_col / 100.0,
			   is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] function_name; 
    delete[] parameter_list; 
    print_html_text_block( *current_stream, T);
    three_cols_html_end( *current_stream, 
			 exp_size > table_width * table_second_col / 100.0,
			 is_empty_comment);
}

void format_variable( const char* signature, 
		      const Text& T, 
		      bool is_typedef = false) {
    char* return_value;
    char* scope;
    char* variable_name; 
    char* rest;    // possibly including assignment

    split_variable_declaration( signature, 
				return_value,
				scope,
				variable_name,
				rest);
    char* formatted_var = convert_ascii_to_html( variable_name);

    double exp_size_ret = 0.0;
    double exp_size     = 0.0;
    if ( return_value) {
        remove_const_ref_pair( return_value);
	exp_size_ret += estimate_html_size( return_value);
    }
 
    three_cols_html_begin( *current_stream, 
			   exp_size_ret > table_width*table_first_col/100.0);
   
    if ( class_name == NULL) {
	// generate a substitution rule for hyperlinking
	*anchor_stream << "[a-zA-Z0-9_]\"" << formatted_var
		       << "\"    { ECHO; }" << endl;
	*anchor_stream << "\"" << formatted_var
		       << "\"[a-zA-Z0-9_]    { ECHO; }" << endl;
	*anchor_stream << '"' << formatted_var
		       << "\"    { fputs( \"<A HREF=\\\""
		       << current_filename 
		       << (is_typedef ? "#Typedef_" :  "#Var_" )
		       << variable_name 
		       << "\\\">" << formatted_var << "</A>\", stdout); }" 
		       << endl;

	// index
	*index_stream << (is_typedef ? sort_key_typedef :
			  sort_key_variable) << '1'; 
	filter_for_index_comment( *index_stream, variable_name);
	*index_stream << "!><UL><LI><VAR>" << formatted_var
		      << "</VAR></UL>" << endl;
    } else {
	// index
	char* scrambled_var = convert_ascii_to_scrambled_html( variable_name);
	*index_stream << (is_typedef ? sort_key_typedef :
			  sort_key_variable) << '1'; 
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, variable_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename 
		      << (is_typedef ? "#Typedef_" :  "#Var_" )
		      << variable_name
		      << "\"><VAR>" << scrambled_template_class_name
		      << scrambled_var << "</VAR></A></UL>" << endl;
	delete[] scrambled_var;
    }
    *current_stream << "<A NAME=\""
		    << (is_typedef ? "Typedef_" :  "Var_" )
		    << variable_name << "\"></A>" << endl;
    // end index

    if ( return_value)
        print_ascii_to_html_spc( *current_stream, return_value);

    // handle function body or operation signature
    // first, estimate size
    if ( scope)
        exp_size += estimate_html_size( scope);
    exp_size += estimate_html_size( variable_name);
    if ( rest)
        exp_size += estimate_html_size( rest);

    // then, do the printing
    three_cols_html_second(
	    *current_stream,
	    exp_size_ret > table_width * table_first_col / 100.0,
	    exp_size > table_width * table_second_col / 100.0
	);
    if ( scope)
        print_ascii_to_html_spc( *current_stream, scope);
    *current_stream << formatted_var;
    if ( rest) {
        *current_stream << ' ';
        print_ascii_to_html_spc( *current_stream, rest);
    }
    *current_stream << ';';

    bool is_empty_comment = is_text_block_empty( T);
    three_cols_html_third( *current_stream, 
			   exp_size > table_width * table_second_col / 100.0,
			   is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] formatted_var; 
    delete[] variable_name; 
    delete[] rest; 
    print_html_text_block( *current_stream, T);
    three_cols_html_end( *current_stream, 
			 exp_size > table_width * table_second_col / 100.0,
			 is_empty_comment);
}

void format_constructor( const char* signature, const Text& T) {
    char* return_value;
    char* scope;
    char* function_name; 
    char* parameter_list; 
    const char* rest;

    split_function_declaration( signature, 
				return_value,
				scope,
				function_name,
				parameter_list,
				rest);

    two_cols_html_begin( *current_stream);
    // first, estimate size
    double exp_size = 0.0;
    if ( scope)
	exp_size += estimate_html_size( scope);
    exp_size += estimate_html_size( creationvariable) 
              + estimate_html_size( template_class_name) 
              + width_per_character;
    int n = 0;
    if ( parameter_list) {
        exp_size += 2 * width_per_character;
        n = separate_parameter_list( parameter_list);
	char* p = parameter_list;
	int m = n;
	while ( m--) {
	    remove_const_ref_pair( p);
	    remove_own_classname( p, template_class_name);
	    exp_size += estimate_html_size( p) + width_per_character;
	    p += strlen( p) + 1;  // skip to next parameter
	}
    }
    exp_size *= stretch_factor;
    if ( exp_size > table_width && parameter_list)
        *current_stream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	                   "<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
	                << indNewline;

    if ( scope)
        print_ascii_to_html_spc( *current_stream, scope);
    print_ascii_to_html_spc( *current_stream, template_class_name);
    *current_stream << ' ';
    print_ascii_to_html_spc( *current_stream, creationvariable);

    if ( parameter_list) {
        *current_stream << " ( ";
	if ( exp_size > table_width)
	    *current_stream << "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
			    << indNewline;
	char* p = parameter_list;
	while ( n--) {
	    print_ascii_to_html_spc( *current_stream, p);
	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_stream << ", ";
		if ( exp_size > table_width)
		    *current_stream << "<BR>" << indNewline;
	    }
	}
        *current_stream << ");";
	if ( exp_size > table_width)
	    *current_stream << "</TD></TR></TABLE>" << indNewline;
    } else
        *current_stream << ';';

    bool is_empty_comment = is_text_block_empty( T);
    two_cols_html_second( *current_stream, is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] function_name; 
    delete[] parameter_list; 
    print_html_text_block( *current_stream, T);
    *current_stream << indNewline << "<P>";
    two_cols_html_end( *current_stream, is_empty_comment);
}

void format_nested_type( const char* nested_type_name, const Text& T) {
    char* formatted_type = convert_ascii_to_html( nested_type_name);
    char* scrambled_type = convert_ascii_to_scrambled_html( nested_type_name);
    two_cols_html_begin( *current_stream);

    // index
    *index_stream << sort_key_nested_type << '1';
    filter_for_index_comment( *index_stream, template_class_name);
    *index_stream << "::";
    filter_for_index_comment( *index_stream, nested_type_name);
    *index_stream << "!><UL><LI><A HREF=\""
		  << current_filename << "#Nested_type_" << nested_type_name
                  << "\"><VAR>" << scrambled_template_class_name
		  << scrambled_type << "</VAR></A></UL>" << endl;
    *current_stream << "<A NAME=\"Nested_type_" << nested_type_name 
		    << "\"></A>" << endl;
    // end index

    // first, estimate size
    double exp_size =   estimate_html_size( template_class_name)
                      + estimate_html_size( nested_type_name) 
                      + 2 * width_per_character;
    exp_size *= stretch_factor;

    print_ascii_to_html_spc( *current_stream, template_class_name);
    *current_stream << "::";
    print_ascii_to_html_spc( *current_stream, nested_type_name);

    bool is_empty_comment = is_text_block_empty( T);
    two_cols_html_second( *current_stream, is_empty_comment);
    print_html_text_block( *current_stream, T);
    *current_stream << indNewline << "<P>";
    two_cols_html_end( *current_stream, is_empty_comment);
    delete[] formatted_type;
    delete[] scrambled_type;
}

void format_enum( const char* signature, const Text& T) {
    char* return_value;
    char* scope;
    char* enum_name; 
    char* parameter_list; 
    const char* rest;

    split_function_declaration( signature, 
				return_value,
				scope,
				enum_name,
				parameter_list,
				rest,
				true);

    char* formatted_enum = convert_ascii_to_html( enum_name);

    two_cols_html_begin( *current_stream);

    if ( class_name == NULL) {
	// generate a substitution rule for hyperlinking
	*anchor_stream << "[a-zA-Z0-9_]\"" << formatted_enum
		       << "\"    { ECHO; }" << endl;
	*anchor_stream << "\"" << formatted_enum
		       << "\"[a-zA-Z0-9_]    { ECHO; }" << endl;
	*anchor_stream << '"' << formatted_enum
		       << "\"    { fputs( \"<A HREF=\\\""
		       << current_filename << "#Enum_" << enum_name << "\\\">"
		       << formatted_enum << "</A>\", stdout); }" 
		       << endl;

	// index
	*index_stream << sort_key_enum << '1';
	filter_for_index_comment( *index_stream, enum_name);
	*index_stream << "!><UL><LI><VAR>" << formatted_enum
		      << "</VAR></UL>" << endl;
    } else {
	// index
	char* scrambled_enum = convert_ascii_to_scrambled_html( enum_name);
	*index_stream << sort_key_enum << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, enum_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Enum_" << enum_name
		      << "\"><VAR>" << scrambled_template_class_name
		      << scrambled_enum << "</VAR></A></UL>" << endl;
	delete[] scrambled_enum;
    }
    *current_stream << "<A NAME=\"Enum_" << enum_name << "\"></A>" << endl;
    // end index

    // first, estimate size
    double exp_size = 0.0;
    if ( scope)
	exp_size += estimate_html_size( scope);
    if ( return_value)
	exp_size += estimate_html_size( return_value);


    exp_size += estimate_html_size( formatted_enum) 
              + 3.0 * width_per_character;
    int n = 0;
    if ( parameter_list) {
        exp_size += 2 * width_per_character;
        n = separate_parameter_list( parameter_list);
	char* p = parameter_list;
	int m = n;
	while ( m--) {
	    exp_size += estimate_html_size( p) + width_per_character;
	    p += strlen( p) + 1;  // skip to next parameter
	}
    }
    exp_size *= stretch_factor;
    if ( exp_size > table_width && parameter_list)
        *current_stream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	                   "<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
	                << indNewline;

    if ( return_value) {
        print_ascii_to_html_spc( *current_stream, return_value);
	*current_stream << ' ';
    }
    if ( scope)
        print_ascii_to_html_spc( *current_stream, scope);
    *current_stream << formatted_enum;
    if ( parameter_list) {
        *current_stream << " { ";
	if ( exp_size > table_width)
	    *current_stream << "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
			    << indNewline;
	char* p = parameter_list;
	while ( n--) {
	    while ( *p && *p <= ' ')
	        ++p;
	    print_ascii_to_html_spc( *current_stream, p);

	    if ( class_name == NULL) {
		// index: print enum tags with their (possible) initializers
		*index_stream << sort_key_enum_tags  << '1' << p 
			      << "!><UL><LI><VAR>";
		print_ascii_to_html_spc( *index_stream, p);
		*index_stream << "</VAR></UL>" << endl;

		// generate a substitution rule for hyperlinking
		// Here, the initializer has to be suppressed
		char* q = p;
		while( *q && *q != '=')
		    ++q;
		while ( q>p && q[-1] <= ' ')
		    --q;
		char c_tmp = *q;
		*q = '\0';
		char *tmp_param = convert_ascii_to_html( p);
		*anchor_stream << "[a-zA-Z0-9_]\"" << tmp_param
			       << "\"    { ECHO; }" << endl;
		*anchor_stream << "\"" << tmp_param
			       << "\"[a-zA-Z0-9_]    { ECHO; }" << endl;
		*anchor_stream << '"' << tmp_param
			       << "\"    { fputs( \"<A HREF=\\\""
			       << current_filename << "#Enum_" 
			       << enum_name << "\\\">"
			       << tmp_param << "</A>\", stdout); }" 
			       << endl;
		delete[] tmp_param;
		*q = c_tmp; // restore initializer
	    } else {
		// index: print enum tags with their (possible) initializers
	        char* scrambled_tag = convert_ascii_to_scrambled_html( p);
		*index_stream << sort_key_enum_tags  << '1' << p 
			      << "!><UL><LI><VAR>";
		print_ascii_to_html_spc( *index_stream, p);
		*index_stream << "</VAR></UL>" << endl;

		*index_stream << sort_key_enum_tags << '1';
		filter_for_index_comment( *index_stream, template_class_name);
		*index_stream << "::";
		filter_for_index_comment( *index_stream, p);
		*index_stream << "!><UL><LI><A HREF=\""
			      << current_filename << "#Enum_" << enum_name
			      << "\"><VAR>" << scrambled_template_class_name
			      << scrambled_tag << "</VAR></A></UL>" << endl;
		delete[] scrambled_tag;
	    }

	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_stream << ", ";
		if ( exp_size > table_width)
		    *current_stream << "<BR>" << indNewline;
	    }
	}
        *current_stream << "};";
	if ( exp_size > table_width)
	    *current_stream << "</TD></TR></TABLE>" << indNewline;
    } else
        *current_stream << ';';

    bool is_empty_comment = is_text_block_empty( T);
    two_cols_html_second( *current_stream, is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] formatted_enum;
    delete[] enum_name;
    delete[] parameter_list; 
    print_html_text_block( *current_stream, T);
    *current_stream << indNewline << "<P>";
    two_cols_html_end( *current_stream, is_empty_comment);
}


void format_struct( const char* signature, const Text& T) {
    char* return_value;
    char* scope;
    char* struct_name; 
    char* parameter_list; 
    const char* rest;

    split_function_declaration( signature, 
				return_value,
				scope,
				struct_name,
				parameter_list,
				rest,
				true);

    char* formatted_struct = convert_ascii_to_html( struct_name);
    two_cols_html_begin( *current_stream);

    if ( class_name == NULL) {
	// generate a substitution rule for hyperlinking
	*anchor_stream << "[a-zA-Z0-9_]\"" << formatted_struct
		       << "\"    { ECHO; }" << endl;
	*anchor_stream << "\"" << formatted_struct
		       << "\"[a-zA-Z0-9_]    { ECHO; }" << endl;
	*anchor_stream << '"' << formatted_struct
		       << "\"    { fputs( \"<A HREF=\\\""
		       << current_filename << "#Struct_" << struct_name 
		       << "\\\">" << formatted_struct << "</A>\", stdout); }" 
		       << endl;
	// index
	*index_stream << sort_key_struct << '1';
	filter_for_index_comment( *index_stream, struct_name);
	*index_stream << "!><UL><LI><VAR>" << formatted_struct
		      << "</VAR></UL>" << endl;
    } else {
	// index
	char* scrambled_struct = convert_ascii_to_scrambled_html( struct_name);
	*index_stream << sort_key_struct << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, struct_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Struct_" << struct_name
		      << "\"><VAR>" << scrambled_template_class_name
		      << scrambled_struct << "</VAR></A></UL>" << endl;
	delete[] scrambled_struct;
    }
    *current_stream << "<A NAME=\"Struct_" << struct_name << "\"></A>" << endl;
    // end index

    // first, estimate size
    double exp_size = 0.0;
    if ( scope)
	exp_size += estimate_html_size( scope);
    if ( return_value)
	exp_size += estimate_html_size( return_value);


    exp_size += estimate_html_size( formatted_struct) 
              + 3.0 * width_per_character;
    int n = 0;
    if ( parameter_list) {
        exp_size += 2 * width_per_character;
        n = separate_parameter_list( parameter_list);
	char* p = parameter_list;
	int m = n;
	while ( m--) {
	    exp_size += estimate_html_size( p) + width_per_character;
	    p += strlen( p) + 1;  // skip to next parameter
	}
    }
    exp_size *= stretch_factor;
    if ( exp_size > table_width && parameter_list)
        *current_stream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	                   "<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
	                << indNewline;

    if ( return_value) {
        print_ascii_to_html_spc( *current_stream, return_value);
	*current_stream << ' ';
    }
    if ( scope)
        print_ascii_to_html_spc( *current_stream, scope);
    *current_stream << formatted_struct;
    if ( parameter_list) {
        *current_stream << " { ";
	if ( exp_size > table_width)
	    *current_stream << "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>" 
			    << indNewline;
	char* p = parameter_list;
	while ( n--) {
	    while ( *p && *p <= ' ')
	        ++p;
	    print_ascii_to_html_spc( *current_stream, p);

	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_stream << ", ";
		if ( exp_size > table_width)
		    *current_stream << "<BR>" << indNewline;
	    }
	}
        *current_stream << "};";
	if ( exp_size > table_width)
	    *current_stream << "</TD></TR></TABLE>" << indNewline;
    } else
        *current_stream << ';';

    bool is_empty_comment = is_text_block_empty( T);
    two_cols_html_second( *current_stream, is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] formatted_struct;
    delete[] struct_name;
    delete[] parameter_list; 
    print_html_text_block( *current_stream, T);
    *current_stream << indNewline << "<P>";
    two_cols_html_end( *current_stream, is_empty_comment);
}


/* Taylored semantic functions used in syntax.y */
/* ============================================ */

void handleMainComment( const Text& T) {
    //  *current_stream << indNewline << "<P>" << endl << indNewline;
    if ( print_html_text_block( *current_stream, T, true))
        *current_stream << indNewline << "<P>" << endl;
}

// void handleComment( const Text& T) {
//     print_html_text_block( *current_stream, T);
//     three_cols_html_end( *current_stream);
// }

// void handleConstructorComment( const Text& T) {
//     print_html_text_block( *current_stream, T);
//     *current_stream << indNewline << "<P>";
//     two_cols_html_end( *current_stream);
// }

void handleChapter(  const Text& T) {
    char* tmp_name = replaceSuffix( in_filename, html_suffix);
    char* new_main_filename = addPrefix( chapter_prefix, tmp_name);
    delete[] tmp_name;
    if ( strcmp( new_main_filename, main_filename) == 0) {
        printErrorMessage( ChapterStructureError);
	delete[] new_main_filename;
	return;
    }
    if ( class_stream != 0) {
        printFootnotes( *class_stream);
        if ( chapter_title)
	    // navigation footer
	    *class_stream << "<HR> Return to chapter: <A HREF=\"" 
			  << main_filename 
			  << "\">" << chapter_title << "</A>" << endl;
        close_html( *class_stream);
	assert_file_write( *class_stream, class_filename);
	delete   class_stream;
	delete[] class_filename;
	class_stream = 0;
	class_filename = 0;
    }
    delete[] chapter_title;
    chapter_title = text_block_to_string( T);
    footnotes        = &main_footnotes;
    footnote_counter = &main_footnote_counter;
    if ( main_stream != &cout) {
        printFootnotes( *main_stream);
        // navigation footer
        *main_stream << "<HR> Next chapter: <A HREF=\"" 
		     << new_main_filename 
		     << "\">" << chapter_title << "</A>" << endl;
        close_html( *main_stream);
	assert_file_write( *main_stream, main_filename);
	delete   main_stream;
	delete[] main_filename;
	main_stream = 0;
	main_filename = 0;
    }
    main_filename = new_main_filename;
    main_stream = open_file_with_path_for_write( tmp_path, main_filename);
    current_stream   = main_stream;
    current_filename = main_filename;
    open_html( *main_stream);

    *main_stream << "<H1>" << chapter_title << "</H1>" << endl;

    // table of contents
    *contents_stream << "    <LI> <A HREF=\"" << main_filename 
		     << "\">" << chapter_title << "</A>" << endl;
}

void handleBiblio(  const Text& T) {
    first_bibitem = true;
    ostream* out = open_file_with_path_for_write( tmp_path, bib_filename);
    istream* in  = open_config_file( cc_biblio_header);
    filter_config_file( *in, *out);

    print_html_text_block( *out, T);
    *out << "</TD></TR>" << endl;

    delete in;
    in = open_config_file( cc_biblio_footer);
    filter_config_file( *in, *out);

    assert_file_write( *out, bib_filename);
    delete in;
    delete out;
}

Buffer* handleCite( const char* l) {
    Buffer* buf = new Buffer;
    buf->add( '[');
    while ( *l != '{' && *l != '[')
        ++l;
    const char* comment = 0;
    const char* end_comment = 0;
    if ( *l == '[') {
        comment = l + 1;
	while ( *l != ']')
	    ++l;
	end_comment = l;
	while ( *l != '{')
	    ++l;
    }
    ++l;
    while( *l != '}') {
        const char* p = l;
	while ( *p != '}' && *p != ',')
	    ++p;
	buf->add( "<A HREF=\"");
	buf->add( bib_filename);
	buf->add( "#Biblio_");
	const char* q = l;
	while ( q != p)
	    buf->add( *q++);
	buf->add( "\">");
	q = l;
	while ( q != p)
	    buf->add( *q++);
	buf->add( "</A>");
	l = p;
	if ( *l == ',') {
	    ++l;
	    buf->add( ", ");
	}
    }
    if ( comment) 
	while ( comment != end_comment)
	    buf->add( *comment++);
    buf->add( ']');
    return buf;
}

// for an empty item name use the key name as item name
Buffer* handleBibItem( const char* key_name, const char* item_name) {
    const char* name = key_name;
    if ( item_name)
        name = item_name;
    Buffer* buf = new Buffer;
    if ( ! first_bibitem) {
        first_bibitem = false;
	buf->add( "</TD></TR>");
    }
    buf->add( "<TR><TD ALIGN=LEFT "
	     "VALIGN=TOP NOWRAP><A NAME=\"Biblio_");
    buf->add( key_name);
    buf->add( "\"></A><B>[");
    buf->add( name);
    buf->add( "]</B></TD><TD ALIGN=LEFT "
	     "VALIGN=TOP>");
    if ( item_name) {
        *anchor_stream << "\"#Biblio_" << key_name
		       << "\"[\"]\">"  << key_name
		       << "\"        { fputs( \"#Biblio_"
		       << key_name << "\\\">" << item_name 
		       << "\", stdout); }" 
		       << endl;
    }
    return buf;
}


void handleSection(  const Text& T) {
    static int section_counter = 1;
    char* section = text_block_to_string( T);
    *current_stream << "<A NAME=\"Section_" << section_counter << "\"></A>"
		    << endl;
    *current_stream << "<H2>" << section << "</H2>" << endl;

    // table of contents
    *contents_stream << "        <UL><LI><A HREF=\"" << current_filename
		     << "#Section_" << section_counter << "\">"
		     << section << "</A></UL>" << endl;

    ++ section_counter;
    delete[] section;
}

void handleLabel( const char* l) {
    /* The lexical processing has removed the parantheses around */
    /* \ref{...} macros from TeX, so here is the correct pattern match */
    /* to find them in the pre-HTML text */
    *anchor_stream << "[\\\\](page)?ref[ \\t]*\"" << l
		   << "\"    { fputs( \"<A HREF=\\\""
		   << current_filename 
		   << "#" << l << "\\\">here</A>\", stdout); }" << endl;    
    // There are two special ref commands defined within the manual
    *anchor_stream << "[\\\\]Chapter[ \\t]*\"" << l
		   << "\"    { fputs( \"<A HREF=\\\""
		   << current_filename 
		   << "#" << l << "\\\">This Chapter</A>\", stdout); }" 
		   << endl;    
    *anchor_stream << "[\\\\]Section[ \\t]*\"" << l
		   << "\"    { fputs( \"<A HREF=\\\""
		   << current_filename 
		   << "#" << l << "\\\">This Section</A>\", stdout); }" 
		   << endl;    
}

void handleText( const Text& T, bool check_nlnl) {
    if ( ! print_html_text_block( *current_stream, T)) {
        if (  check_nlnl) {
	    int count = 0;
	    InListFIter< TextToken> words( (Text&)T);
	    ForAll( words) {
		if ( (*words).string[0] == '\n')
		    ++count;
	    }
	    if ( count > 1)
		*current_stream << "<P>\n" << endl;
	    else if ( count > 0)
	        *current_stream << '\n';
	   else
	        *current_stream << ' ';
	} else
	    *current_stream << ' ';
    }
}

void handleBuffer( const Buffer& B) {
    *current_stream << B.string();
}

void handleTextToken( const TextToken& TT) {
    *current_stream << TT.string;
}

void handleString( const char* s) {
    *current_stream << s;
}

void handleChar( char c) {
    *current_stream << c;
}


void handleClasses( const char* classname, const char* template_cls) {
    if ( template_cls)
        template_class_name  = newstr( template_cls);
    else
        template_class_name  = newstr( classname);
    assert( template_class_name);
    char* t_tmp_name = convert_ascii_to_html( template_class_name);
    formatted_template_class_name = new char[ strlen( t_tmp_name) + 12];
    strcpy( formatted_template_class_name, "<VAR>");
    strcat( formatted_template_class_name, t_tmp_name);
    strcat( formatted_template_class_name, "</VAR>");
    delete[] t_tmp_name;
    t_tmp_name = convert_ascii_to_scrambled_html( template_class_name);
    scrambled_template_class_name = new char[ strlen( t_tmp_name) + 12];
    strcpy( scrambled_template_class_name, "<VAR>");
    strcat( scrambled_template_class_name, t_tmp_name);
    strcat( scrambled_template_class_name, "</VAR>");
    delete[] t_tmp_name;

    if ( class_stream != 0) {
        printFootnotes( *class_stream);
        // navigation footer
        *class_stream << "<HR> Next class declaration: "
		      << formatted_template_class_name << endl;
        close_html( *class_stream);
	assert_file_write( *class_stream, class_filename);
	delete   class_stream;
	delete[] class_filename;
	class_stream = 0;
    }
    class_name           = newstr( classname);
    footnotes        = &class_footnotes;
    footnote_counter = &class_footnote_counter;

    char *tmp_name = convert_ascii_to_html( classname);
    formatted_class_name = new char[ strlen( tmp_name) + 12];
    strcpy( formatted_class_name, "<VAR>");
    strcat( formatted_class_name, tmp_name);
    strcat( formatted_class_name, "</VAR>");

    class_filename = addSuffix( classname, html_suffix);
    class_stream = open_file_with_path_for_write( tmp_path, class_filename);
    open_html( *class_stream);
    // Make a hyperlink in the chapter to the class file.
    if ( main_stream != &cout) {
        *main_stream  << "<UL><LI>\nClass declaration for "
		      << formatted_template_class_name
		      << ".</UL>\n" << endl;
    }

    // table of contents
    *contents_stream << "        <UL><LI> Class declaration of " 
		     << formatted_template_class_name << "</UL>" << endl;

    // index
    *index_stream << sort_key_class << '1';
    filter_for_index_comment( *index_stream, template_class_name);
    *index_stream << "!><UL><LI>" << formatted_template_class_name 
		  << "</UL>" << endl;

    current_stream   = class_stream;
    current_filename = class_filename;

    // generate a substitution rule for hyperlinking
    *anchor_stream << '"' << class_filename
		   << "\"    { ECHO; }" << endl;
    *anchor_stream << "[a-zA-Z0-9_]\"" << tmp_name
		   << "\"    { ECHO; }" << endl;
    *anchor_stream << "\"" << tmp_name
		   << "\"[a-zA-Z0-9_]    { ECHO; }" << endl;
    if ( template_cls) {
        *anchor_stream << '"' << tmp_name
		       << "\"[ ]*\"&lt;\"    {\n"
		       << "        fputs( \"<A HREF=\\\""
		       << class_filename << "\\\">\", stdout);\n" 
		       << "        nesting = 1;\n"
		       << "        yymore();\n"
		       << "        BEGIN( PARAMMODE); }\n"
		       << endl;
    }
    *anchor_stream << '"' << tmp_name
		   << "\"    { fputs( \"<A HREF=\\\""
		   << class_filename << "\\\">"
		   << tmp_name << "</A>\", stdout); }" 
		   << endl;
    delete[] tmp_name;
    creationvariable = newstr( "this");
    formatted_creationvariable = newstr( "<VAR>this</VAR>");
}

void handleClass( const char* classname) {
    handleClasses( classname, 0);
}

void handleClassEnd( void) {
    delete[] class_name;
    delete[] template_class_name;
    delete[] formatted_class_name;
    delete[] formatted_template_class_name;
    delete[] scrambled_template_class_name;
    class_name = 0;
    formatted_class_name = 0;
    template_class_name = 0;
    formatted_template_class_name = 0;
    scrambled_template_class_name = 0;
    /* ...  Hack to implement the link from one class to the next class
    close_html( *class_stream);
    if ( ! *class_stream) {
        cerr << ' ' << endl 
	     << prog_name << ": error: cannot write to file `" 
	     <<  class_filename << "'." << endl;
	exit(1);
    }
    delete   class_stream;
    delete[] class_filename;
    class_stream = 0;
    class_filename = 0;
    ... */
    current_stream   = main_stream;
    current_filename = main_filename;
}

void handleClassTemplate( const char* classname) {
    char* s = (char *)classname;
    while ( *s != 0 && *s != '<') s++;
    if ( *s == 0)
        printErrorMessage( TemplateParamExpectedError);
    char c_tmp = *s;
    *s = 0;
    char *classname_tmp = newstr( classname);
    *s = c_tmp;
    handleClasses( classname_tmp, classname);
    delete[] classname_tmp;
}

void handleClassTemplateEnd( void) {
    handleClassEnd();
}


void handleDeclaration( const char* ) {}

void handleMethodDeclaration( const char* decl, const Text& T) {
    decl = handle_template_layout( *current_stream, decl, true);
    format_function( true, decl, T);
}

void handleFunctionDeclaration( const char* decl, const Text& T) {
    decl = handle_template_layout( *current_stream, decl, true);
    format_function( false, decl, T);
}

void handleConstructorDeclaration( const char* decl, const Text& T) {
    decl = handle_template_layout( *current_stream, decl, false);
    format_constructor( decl, T);
}

void handleFunctionTemplateDeclaration( const char* ,
					const char* decl,
					const Text& T) {
    handleFunctionDeclaration( decl, T);
}

void handleVariableDeclaration( const char* decl, const Text& T) {
    format_variable( decl, T);
}

void handleTypedefDeclaration( const char* decl, const Text& T) {
    format_variable( decl, T, true);
}

void handleNestedTypeDeclaration( const char* decl, const Text& T) {
    format_nested_type( decl, T);
}

void handleEnumDeclaration( const char* decl, const Text& T) {
    format_enum( decl, T);
}

void handleStructDeclaration( const char* decl, const Text& T) {
    decl = handle_template_layout( *current_stream, decl, false);
    format_struct( decl, T);
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
    int i;
    int nParameters = 0;
    char *parameters[ MaxParameters + 1];
 
    Switch help_switch  = NO_SWITCH;
    Switch dummy_switch;
    for (i = 1; i < argc; i++) {

        /* check switches */
        detectSwitch( trace_switch, "trace");
	    yydebug = 1;
        endDetect();
        detectSwitch( line_switch, "line");
        endDetect();

        detectSwitch( dummy_switch, "date");
            i++;
            if ( i < argc)
	        manual_date = argv[i];
	    else {
	        cerr << "error: option -date needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "release");
            i++;
            if ( i < argc)
	        manual_release = argv[i];
	    else {
	        cerr << "error: option -release needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "title");
            i++;
            if ( i < argc)
	        manual_title = argv[i];
	    else {
	        cerr << "error: option -title needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "author");
            i++;
            if ( i < argc)
	        manual_author = argv[i];
	    else {
	        cerr << "error: option -author needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( config_switch, "config");
            i++;
            if ( i < argc) {
	        config_path = argv[i];
		if ( config_path[ strlen( config_path) - 1] != '/') {
		    cerr << "error: option -config: a path must terminate "
		            "with a /" << endl;
		    nParameters = ErrParameters;
		}
	    } else {
	        cerr << "error: option -config needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( dummy_switch, "tmp");
            i++;
            if ( i < argc) {
	        tmp_path = argv[i];
		if ( tmp_path[ strlen( tmp_path) - 1] != '/') {
		    cerr << "error: option -tmp: a path must terminate "
		            "with a /" << endl;
		    nParameters = ErrParameters;
		}
	    } else {
	        cerr << "error: option -tmp needs an additional parameter"
		     << endl;
	        nParameters = ErrParameters;
	    }
        endDetect();
        detectSwitch( warn_switch, "warn");
        endDetect();
        detectSwitch( macro_switch, "macro");
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
	cerr << prog_name << " " << prog_release << " (c) Lutz Kettner" 
	     << endl;
        cerr << "Usage: " << prog_name << " [<options>] <TeX-files...>" 
	     << endl;
        cerr << "       -date    <text>    set a date for the manual." << endl;
        cerr << "       -release <text>    set a release number for the "
                                          "manual." << endl;
        cerr << "       -title   <text>    set a title text for the manual." 
	     << endl;
        cerr << "       -author  <text>    set an author address (email) for "
	                                  "the manual." << endl;
        cerr << "       -config  <dir/>    set the path where to find the "
                                          "config files." << endl;
        cerr << "       -tmp     <dir/>    set the path where to put the "
                                          "output files." << endl;
        cerr << "       -warn              warn about unknown macros" << endl;
        cerr << "       -macro             trace macro definitions" << endl;
        cerr << "       -trace             set `yydebug' for bison to true"
             << endl;
        cerr << "       -line              echo currently parsed line "
                                          "numbers to cerr" << endl;
        cerr << "(TeX-files with suffix:  .tex  or  .bbl  possible.)"  
	     << endl;
        exit(1);
    }


    // Initialization
    tag_defaults();
    main_stream   = &cout;

    // Filter config files for the index.
    copy_config_file( cc_index_header);
    copy_config_file( cc_index_footer);

    // Prepare several streams:
    // anchor_stream = new ofstream( anchor_filename, ios::out | ios::app);
    anchor_stream   = open_file_with_path_for_write( tmp_path, 
						     anchor_filename);
    contents_stream = open_file_with_path_for_write( tmp_path, 
						     contents_filename);
    copy_and_filter_config_file( cc_toc_header, *contents_stream);
    index_stream    = open_file_with_path_for_write( tmp_path, index_filename);
    write_headers_to_index( *index_stream);
    
    for ( i = 0; i < nParameters; i++) {
	FILE* in;
	if ( (in = fopen( parameters[i], "r")) == NULL) {
	    cerr << ' ' << endl 
		 << prog_name << ": error: cannot open infile `" 
		 << parameters[i] << "'." << endl;
	    exit(1);
	}
	in_filename = parameters[i];

	main_stream   = &cout;
	current_stream   = main_stream;
	current_filename = main_filename;

	init_scanner( in);
	yyparse();
	fclose( in);

	if ( class_stream != 0) {
            printFootnotes( *class_stream);
	    close_html( *class_stream);
	    assert_file_write( *class_stream, class_filename);
	    delete   class_stream;
	    delete[] class_filename;
	    class_stream = 0;
	    class_filename = 0;
	}
	assert_file_write( *main_stream, main_filename);
	if ( main_stream != &cout) {
            footnotes        = &main_footnotes;
            footnote_counter = &main_footnote_counter;
            printFootnotes( *main_stream);
	    close_html( *main_stream);
	    assert_file_write( *main_stream, main_filename);
	    delete   main_stream;
	    delete[] main_filename;
	    main_stream = &cout;
	    main_filename = "<cout>";
	}
    }

    assert_file_write( *index_stream, index_filename);
    delete index_stream;

    copy_and_filter_config_file( cc_toc_footer, *contents_stream);
    assert_file_write( *contents_stream, contents_filename);
    delete contents_stream;

    assert_file_write( *anchor_stream, anchor_filename);
    delete anchor_stream;
    delete[] chapter_title;

    cout << endl;
    return 0;
}

// EOF //

