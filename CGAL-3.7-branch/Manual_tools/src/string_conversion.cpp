/**************************************************************************
 
  string_conversion.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : String conversion functions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <string_conversion.h>
#include <stdlib.h>
#include <ctype.h>
#include <output.h>
#include <input.h>
#include <sstream>
#include <config.h>
#include <error.h>
#include <macro_dictionary.h>

// New style conversion routines
// =======================================

template< typename T >
string anytype_to_string( T t ) {
   ostringstream out;
   out << t;
   return out.str();
}

string   int_to_string( int   i ) { return anytype_to_string( i ); }
string float_to_string( float f ) { return anytype_to_string( f ); }

// Roman numbers up to 109
static const char * const roman_numbers[110] = {
    "0", "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix",
    "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii", "xix",
    "xx","xxi","xxii","xxiii","xxiv","xxv","xxvi","xxvii","xxviii","xxix",
    "xxx","xxxi","xxxii","xxxiii","xxxiv","xxxv","xxxvi","xxxvii",
        "xxxviii","xxxix",
    "xl","xli","xlii","xliii","xliv","xlv","xlvi","xlvii","xlviii","xlix",
    "l","li","lii","liii","liv","lv","lvi","lvii","lviii","lix",
    "lx","lxi","lxii","lxiii","lxiv","lxv","lxvi","lxvii","lxviii","lxix",
    "lxx","lxxi","lxxii","lxxiii","lxxiv","lxxv","lxxvi","lxxvii","lxxviii",
        "lxxix",
    "lxxx","lxxxi","lxxxii","lxxxiii","lxxxiv","lxxxv","lxxxvi","lxxxvii",
        "lxxxviii","lxxxix",
    "xc","xci","xcii","xciii","xciv","xcv","xcvi","xcvii","xcviii","xcix",
    "c","ci","cii","ciii","civ","cv","cvi","cvii","cviii","cix"
};


/* An object storing the current font */
/* ================================== */
/* It is only used within CCMode at the moment */

Font current_font = rm_font;

/* HTML Tags to set and reset a specific font. */
const char* html_font_opening[ end_font_array] = {
    "",
    "<TT>",
    "<B>",
    "<I>",
    "<I>",
    "<TT>",
    "<TT>",
    "<VAR>",
    "<span class=\"math\">"
};
const char* html_font_closing[ end_font_array] = {
    "",
    "</TT>",
    "</B>",
    "</I>",
    "</I>",
    "</TT>",
    "</TT>",
    "</VAR>",
    "</span>"
};

const int max_tag_size = 20;
char  font_tag_buffer[max_tag_size];


// Returns the roman digit representation for a number i, 0 <= i <= 109.
string int_to_roman_string( int i) {
    if ( i < 1 || i > 109) {
	printErrorMessage( RomansOutOfBoundsError);
        return string( "[Roman digits out of bounds]");
    }
    return string( roman_numbers[i]);
}

void remove_leading_spaces( string& s) {
    size_t i = 0;
    while ( i < s.size() && ( isspace( s[i]) || s[i] == SEPARATOR))
	++i;
    if ( i > 0)
	s.replace( 0, i, "");
}

void remove_trailing_spaces( string& s) {
    if ( s.empty())
	return;
    int i = s.size() - 1;
    while ( i > 0 && ( isspace( s[i]) || s[i] == SEPARATOR))
	--i;
    if ( i < int(s.size()) - 1)
	s.replace( i+1, s.size() - i - 1, "");
}

void crop_string( string& s) {
    remove_leading_spaces( s);
    remove_trailing_spaces( s);
}

void compress_spaces_in_string( string& s) {
    for ( size_t i = 0; i < s.size(); ++i) {
	if ( isspace( s[i])) {
	    size_t k = 1;
	    while ( i + k < s.size() && isspace( s[ i + k]))
		++k;
	    s.replace( i, k, " ");
	}
    }
}

void remove_separator( string& s) {
    for ( size_t i = 0; i < s.size();) {
	if ( s[i] == SEPARATOR) {
	    size_t k = 1;
	    while ( i + k < s.size() && s[ i + k] == SEPARATOR)
		++k;
	    s.replace( i, k, "");
	} else
	    ++i;
    }
}

// Replaces the <> around any template parameters with -- since < and > cannot
// be used in file names under M$Windows.  Also replaces colons (:)  by -'s
// since this character is also not allowed by M$ and commas and space by -'s 
// since these cause problems as well.
string replace_template_braces_and_colons( string name) {
    for ( size_t i = 0; i < name.size(); ++i) {
	if ( name[i]==':' || 
             name[i]==',' || 
             name[i]==' ' || 
             name[i]=='<' ||
             name[i]=='>' /*||
             name[i]=='(' ||
             name[i]==')' */ ) 
        {
	    name.replace(i,1,"-");
	}
    }
    return name;
}

// Replaces all *'s in name by the string "_star" since the * causes problems
// when creating a file (i.e., a message "$f: Ambiguous" is generated and no
// file is created when applying the anchor filter)
string replace_asterisks( string name) {
   for ( size_t i = 0; i < name.size(); ++i) {
      if ( name[i] == '*')
      {
         name.replace(i,1,"_star");
         i +=4;
      }
   }
   return name;
}

// Removes the quoted font changing commands used in CCMode: |I|, |B| ...
string remove_font_commands( string name) {
    for ( size_t i = 0; i + 2 < name.size(); ++i) {
	if ( name[i]=='|' && isupper(name[i+1]) && name[i+2]=='|') {
	    name.replace(i,3,"");
	    --i;
	}
    }
    return name;
}

string remove_suffix( string name) {
    string::size_type i = name.rfind( '.');
    if ( i != string::npos)
	name.replace( i, name.size() - i, "");
    return name;
}

string basename_string( string name) {
    string::size_type i = name.rfind( '/');
    if ( i != string::npos)
	name.replace( 0, i+1, "");
    return name;
}

string rootname_string( string name) {
    name = remove_suffix(name);
    string::size_type i = name.rfind( '/');
    if ( i != string::npos)
	name.replace( 0, i+1, "");
    return name;
}

string path_string( string name) { // either empty or includes a trailing '/'
    // remove leading './' component
    if ( name[0] == '.' && name[1] == '/')
        name.replace(0,2,"");
    string::size_type i = name.rfind( '/');
    if ( i == string::npos)
	return string();
    name.replace( i+1, name.size() - i - 1, "");
    return name;
}

string uppath_string( string path) { // '../../' type of path reversing 'path'
    std::size_t len = std::count( path.begin(), path.end(), '/');
    // absolute paths imply an empty uppath.
    if( len > 1 && path[0] == '/' )
        return string("");
    string result( "");
    for ( std::size_t i = 0; i < len; ++i) {
        result += "../";
    }
    return result;
}

string suffix_string( string name) {
    string::size_type i = name.rfind( '.');
    if ( i != string::npos)
	return name.substr( i+1);
    return string();
}

void assert_trailing_slash_in_path( string& s) {
    if ( s.size() > 0 && s.at( s.size() - 1) != '/') {
	s += '/';
    }
}

// Quoted strings use C string notation with \ as escape symbol.
// Replaced sequences are: \\, \n, \t, \{, \}.
string convert_quoted_string( string s) {
    size_t i = 0;
    while (  i < s.size()) {
	if ( s[i] == SEPARATOR)
	    s.replace( i, 1, "");
	else {
	    if ( s[i] == '\\' && i + 1 < s.size()) {
		switch (s[i+1]) {
		case '\\':
		case '{':
		case '}':
		    s.replace( i, 1, "");
		    break;
		case 'n':
		    s.replace( i, 2, "\n");
		    break;
		case 't':
		    s.replace( i, 2, "\t");
		    break;
		default:
		    break;
		}
	    }
	    ++i;
	}
    }
    return s;
}

// Quoted strings use C string notation with \ as escape symbol.
// Replaced sequences are: \\, \n, \t, \{, \}.
// Makes SEPARATOR Explictly visible for debugging, see \lciDump.
string convert_quoted_string_seps( string s) {
    size_t i = 0;
    while (  i < s.size()) {
	if ( s[i] == SEPARATOR) {
	    s.replace( i, 1, "^A");
            i += 2;
	} else {
	    if ( s[i] == '\\' && i + 1 < s.size()) {
		switch (s[i+1]) {
		case '\\':
		case '{':
		case '}':
		    s.replace( i, 1, "");
		    break;
		case 'n':
		    s.replace( i, 2, "\n");
		    break;
		case 't':
		    s.replace( i, 2, "\t");
		    break;
		default:
		    break;
		}
	    }
	    ++i;
	}
    }
    return s;
}

// Expands '"' and \ symbols with quotes. Replaces \n with \\n, \t with \\t.
string convert_to_C_printable( string s) {
    size_t i = 0;
    while (  i < s.size()) {
	switch (s[i]) {
	case SEPARATOR:
	    s.replace( i, 1, "");
	    break;
	case '\\':
	case '"':
	    s.insert( i, "\\");
	    i += 2;
	    break;
	case '\n':
	    s.replace( i, 1, "\\n");
	    i += 2;
	    break;
	case '\t':
	    s.replace( i, 1, "\\t");
	    i += 2;
	    break;
	default:
	    ++i;
	    break;
	}
    }
    return s;
}

// the following functions are only needed for `wrap_anchor'

/* find the next occurance of c in s or the '\0' character */
const char* find_char( const char* s, char c) {
    while( *s && *s != c )
        ++s;
    return s;
}

/* reverse find of the prev occurance of c in s or the begining position p */
const char* rfind_char( const char* s, const char* p, char c) {
    while( s != p && *s != c )
        --s;
    return s;
}

/* print a HTML tag up to and including the '>' */
void print_tag( const char* s, ostream& out ) {    
    while( *s && *s != '>' ) {
        out << *s;
        ++s;
    }
    if( *s )
        out << *s;
}   

/* wraps an anchor around a body. Checks for font changing tags. */
void wrap_anchor( const string& cc_url, const string& cc_body, ostream& out ) {
    const char* url  =  cc_url.c_str();
    const char* body = cc_body.c_str();
    const char* tag_begin = 0;
    const char* tag_end = 0;
    const char* s = find_char( body, '<');
    if( *s ) {
        ++s;
        if ( *s == '/')
            tag_begin = s + 1;
    }
    s = find_char( body, '\0');
    --s;
    s = rfind_char( s, body, '<');
    if( *s ==  '<' ) {
        ++s;
        if ( *s != '/')
            tag_end = s;
    }
    if( tag_begin ) {
        out << "</";
        print_tag( tag_begin, out );
    }
    out << "<A HREF=\"" << url << "\">";
    if( tag_begin ) {
        out << '<';
        print_tag( tag_begin, out );
    }
    out << body;
    if( tag_end ) {
        out << "</";
        print_tag( tag_end, out );
    }
    out << "</A>";
    if( tag_end ) {
        out << '<';
        print_tag( tag_end, out );
    }
}



// read a file into a string
// ===========================
void append_file_to_string( const string& name, string& s) {
//    istream* in = open_file_for_read( name.c_str());
    istream* in = open_file_for_read_w_input_dirs( name.c_str());
    char c;
    while( in->get(c))
        s += c;
}

void open_file_to_string( const string& name, string& s) {
    istream* in = open_file_for_read( name.c_str());

    char c;
    while( in->get(c))
        s += c;
    delete in;
}



// Old style conversion routines
// =======================================

static bool encode_backslash_flag = true;

void encode_backslash( bool b ) {
  encode_backslash_flag = b;
}

bool
is_html_multi_character( char c ) {
    if( c == SEPARATOR || c == '"' || c == '&' || c == '<' || c == '>' )
      return true;
    return encode_backslash_flag && c == '\\';
}

static char multi_char_default[2];

const char* html_multi_character( char c) {
    switch ( c) {
    case SEPARATOR: return "";
    case '"': return "&quot;";
    case '&': return "&amp;";
    case '<': return "&lt;";
    case '>': return "&gt;";
    case '\\': return "&#92;";
    default:  
	multi_char_default[0] = c;
	multi_char_default[1] = '\0';
    }
    return multi_char_default;
}

void print_ascii_to_html( ostream& out, const char* txt) {
    if (txt == NULL)
	return;
    while( *txt) {
	if ( *txt == '|' && isupper(txt[1]) && txt[2] == '|') {
	    out << new_remember_font( txt[1]);
	    txt += 2;
	} else if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    out << *txt;
	++txt;
    }
}

void print_ascii_len_to_html( ostream& out, const char* txt, int n) {
    if (txt == NULL)
	return;
    while( n) {
	if ( *txt == '|' && isupper(txt[1]) && txt[2] == '|') {
	    out << new_remember_font( txt[1]);
	    txt += 2;
	} else if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    out << *txt;
	++txt;
	n--;
    }
}

// This version eliminates multiple spaces.
void print_ascii_to_html_spc( ostream& out, const char* txt) {
    if (txt == NULL)
	return;
    while( *txt) {
	if ( *txt == '|' && isupper(txt[1]) && txt[2] == '|') {
	    out << new_remember_font( txt[1]);
	    txt += 2;
	} else if ( is_html_multi_character( *txt))
	    out << html_multi_character( *txt);
	else
	    if ( *txt > ' ' || (*txt > '\0' && txt[1] > ' '))
	        out << *txt;
	++txt;
    }
}

int strlen_ascii_to_html( const char* txt) {
    if (txt == NULL)
	return 0;
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


int strlen_fontified_ascii_to_html( const char* txt) {
    if (txt == NULL)
	return 0;
    int len = 0;
    while( *txt) {
	if ( *txt == '|' && isupper(txt[1]) && txt[2] == '|') {
	    len += 13; // upper bound, see </MATH>.
	    txt += 2;
	} else if ( is_html_multi_character( *txt))
	    len += strlen( html_multi_character( *txt));
	else
	    len++;
	++txt;
    }
    return len;
}

char* convert_fontified_ascii_to_html( const char* txt) {
    if ( txt == NULL) {
        char *q = new char[1];
	q[0] = '\0';
	return q;
    }
    char* s = new char[ strlen_fontified_ascii_to_html( txt) + 1];
    char* p = s;
    while( *txt) {
	if ( *txt == '|' && isupper(txt[1]) && txt[2] == '|') {
	    const char* q = new_remember_font( txt[1]);
	    while ( *q)
	        *p++ = *q++;
	    txt += 2;
	} else if ( is_html_multi_character( *txt)) {
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

int strlen_for_makeindex( const char* txt) {
    if (txt == NULL)
	return 0;
    int len = 0;
    while( *txt) {
	if ( *txt == '|' || *txt=='!' || *txt=='@' ) {
	    len+=2;   //add "
	    ++txt;
	} else {
            ++len;
	    ++txt;
        }  
    }
    return len;
}

char* convert_indexentry_for_makeindex(const char* txt) {
  if ( txt == NULL) {
        char *q = new char[1];
	q[0] = '\0';
	return q;
    }
    char* s = new char[ strlen_for_makeindex( txt) + 1];
    char* p = s;
    while( *txt) {
	if ( *txt == '|'  || *txt=='!' || *txt=='@') {
            *p++ = '"';
            *p++ = *txt++;
	} else
	    *p++ = *txt++;
    }
    *p = '\0';
    return s;
}

string string_to_lower( string s ) {
    for( string::iterator it = s.begin(); it != s.end(); ++it ) {
        *it = tolower( *it );
    }
    return s;
}

Font latex_font_to_Font( string font ) {
    font = string_to_lower( font );
    if( font == "\\it" )
        return it_font;
    else if( font == "\\bf" )
        return bf_font;
    
    return rm_font;
}
   
string convert_C_to_html( const char* txt) {
    Font ccFont         = latex_font_to_Font( macroX( "\\ccFont" ) );
    string open_ccFont  = html_font_opening[ ccFont ];
    string close_ccFont = html_font_closing[ ccFont ];
    current_font = ccFont;
    char* tmp = convert_fontified_ascii_to_html( txt );
    string retval = open_ccFont + string(tmp) + close_ccFont;
    delete[] tmp;
    return retval;
}



/* Filter decl. before writing to the index comment */
/* ================================================ */

void filter_for_index_comment( ostream& out, const char* text) {
    if ( text == 0)
	return;
    while( *text) {
        switch (*text) {
	case '|':
	    if ( isupper(text[1]) && text[2] == '|')
		text+=2;
	    else
		out << *text;
	    break;
	case '<':
	    out << '(';	  
	    break;
	case '>':
	    out << ')';	  
	    break;
	case '&':
	    out << '_';	  
	    break;
	case ' ':
	case '\n':
	case '\t':
	case '\r':
	    break;
	default:
	    out << *text;
	}
        text++;
    }
}

string filter_for_index_comment( string s) {
    for ( size_t i = 0; i < s.size(); ++i) {
        switch ( s[i]) {
	case '|':
	    if ((i < s.size() + 2) && isupper(s[i+1]) && s[i+2] == '|')
		s.replace( i--, 3, "");
	    break;
	case '<':
	    s.replace( i, 1, "(");
	    break;
	case '>':
	    s.replace( i, 1, ")");
	    break;
	case '&':
	    s.replace( i, 1, "_");
	    break;
	case ' ':
	case '\n':
	case '\t':
	case '\r':
	    s.replace( i, 1, "");
	    --i;
	    break;
	}
    }
    return s;
}

// Filter characters that might not be allowed in hyperlink anchors.
void filter_for_index_anchor( ostream& out, const char* text) {
    if ( text == 0)
	return;
    while( *text) {
        switch ( *text) {
	case '|':
	    if ( isupper(text[1]) && text[2] == '|')
		text+=2;
	    else
		out << '_';
	    break;
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
	case '\n':
	case '\t':
	case '\r':
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



const char* font_changing_tags( Font old_font, Font new_font) {
    CC_Assert( old_font > unknown_font && old_font < end_font_array);
    CC_Assert( new_font > unknown_font && new_font < end_font_array);
    if ( old_font == new_font)
	*font_tag_buffer = '\0';
    else {
	strcpy( font_tag_buffer, html_font_closing[ old_font]);
	strcat( font_tag_buffer, html_font_opening[ new_font]);
    }
    return font_tag_buffer;
}

const char* new_font_tags( Font new_font) {
    const char* s = font_changing_tags( current_font, new_font);
    current_font = new_font;
    return s;
}

const char* lazy_new_font_tags( Font new_font) {
    char* s = font_tag_buffer;
    *s = '\0';
    if ( current_font == unknown_font && new_font == unknown_font)
	return s;
    if ( current_font == unknown_font)
	strcpy( s, html_font_opening[ new_font]);
    else if ( new_font == unknown_font)
	strcpy( s, html_font_closing[ current_font]);
    else
	return new_font_tags( new_font);
    current_font = new_font;
    return s;
}


/* The following is used to remember the font that was used in the */
/* previous column and to restore it in the next column. */
Font remember_font = it_font;  // default for program code.

const char* new_remember_font( char c) {
    switch (c) {
    case 'T':
	return new_remember_font( tt_font);
    case 'I':
	return new_remember_font( it_font);
    case 'B':
	return new_remember_font( bf_font);
    default:
	font_tag_buffer[0] = '|';
	font_tag_buffer[1] = c;
	font_tag_buffer[2] = '|';
	font_tag_buffer[3] = '\0';
    }
    return font_tag_buffer;
}


// EOF //

