/**************************************************************************
 
  cpp_formatting.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : C++ Formatting functions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <cpp_formatting.h>
#include <string_conversion.h>
#include <macro_dictionary.h>
#include <internal_macros.h>
#include <html_config.h>
#include <html_error.h>
#include <stdlib.h>
#include <ctype.h>
#include <output.h>

// New style conversion routines
// =======================================

string class_name;
string template_class_name;

// Old style conversion routines
// =======================================


/* table size and font size constants */
/* ================================== */
const int table_width             = 650; // absolute
const int table_long_param_indent = 50;  // absolute
const int table_first_col         = 25;  // in percent
const int table_second_col        = 25;  // in percent
const int table_third_col         = 50;  // in percent
const int table_2c_first_col      = 30;  // in percent
const int table_2c_second_col     = 70;  // in percent

const double width_per_character  = 5.5;

// This factor is multiplied to the actual width of an C++ declaration
// right before the test for multiple lines formatting occurs.
// A greater value forces declarations to be printed in multiple lines.
double stretch_factor = 1.6;


/* HTML generating functions */
/* ========================= */

double estimate_html_size( const char* s) {
    if ( s == 0)
	return 0;
    int n = 0;
    while ( *s) {
	if ( *s == '|' && isupper(s[1]) && s[2] == '|') {
	    s += 2;
	} else  if ( *s > ' ' || (*s > '\0' && s[1] > ' '))
	    n++;
        ++s;
    }
    return width_per_character * n; 
}

inline double estimate_html_size( const string& s) {
    return estimate_html_size( s.c_str());
}

void three_cols_html_begin( ostream& out, bool big_col1) {
    current_font = it_font;
    out << indent << indent << ind_newline
	<< "<!3><TABLE BORDER=0 CELLSPACING=2 CELLPADDING=0 WIDTH="
	<< table_width << ">" << ind_newline 
	<< "<TR><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_first_col + (big_col1 ? (table_second_col+table_third_col) :0)
	<< "%" << ( big_col1 ? " COLSPAN=3>" : " NOWRAP>")
	<< ind_newline << "<I><NOBR>" << outdent << ind_newline;
}

void three_cols_html_premature_end( ostream& out) {
    out << indent << ind_newline << store_remember_font() << "</I></NOBR>" 
	<< ind_newline << "</TD></TR>" << ind_newline 
	<< "</TABLE><!3>" << outdent << outdent << ind_newline;
}

void three_cols_html_second( ostream& out, bool big_col1, bool big_col2) {
    out << indent << ind_newline << store_remember_font() << "</I></NOBR>" 
	<< ind_newline << "</TD>";
    if ( big_col1)
        out << "</TR><TR><TD WIDTH=" << table_first_col << "% NOWRAP></TD>";
    out << "<TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_second_col + ( big_col2 ? table_third_col : 0)
	<< "% NOWRAP" << ( big_col2 ? " COLSPAN=2>" : ">")
	<< ind_newline << "<I><NOBR>" << get_remember_font() << outdent 
	<< ind_newline;
}

void three_cols_html_third( ostream& out, bool big_col2, bool empty_col3) {
    out << indent << ind_newline << store_remember_font() << "</I></NOBR>" 
	<< ind_newline;
    if ( ! big_col2)
	out << "</TD>";
    if ( ! empty_col3) {
	if ( big_col2)
	    out << "</TR><TR><TD WIDTH=" << table_first_col 
	        << "% NOWRAP></TD><TD WIDTH=" << table_second_col 
		<< "% NOWRAP></TD>";
	out << "<TD ALIGN=LEFT VALIGN=TOP WIDTH=" << table_third_col << "%>";
    }
    out << outdent << ind_newline;
}

void three_cols_html_end( ostream& out, bool empty_col3) {
    out << indent << ind_newline;
    if ( ! empty_col3)
        out << "</TD>";
    out << "</TR>" << ind_newline 
	<< "</TABLE><!3>" << outdent << outdent << ind_newline;
}

void three_cols_html_new_closing( ostream& out) {
    out << outdent;
}

void two_cols_html_begin( ostream& out) {
    current_font = it_font;
    out << indent << indent << ind_newline
	<< "<!2><TABLE BORDER=0 CELLSPACING=2 CELLPADDING=0 WIDTH="
	<< table_width << ">" << ind_newline 
	<< "<TR><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	<< table_2c_first_col + table_2c_second_col << "% NOWRAP COLSPAN=2>"
	<< ind_newline << "<I><NOBR>" << outdent << ind_newline;
}

void two_cols_html_premature_end( ostream& out) {
    out << indent << ind_newline << store_remember_font() <<"</I></NOBR>" 
	<< ind_newline << "</TD></TR>" << ind_newline 
	<< "</TABLE><!2>" << outdent << outdent << ind_newline;
}

void two_cols_html_second( ostream& out, bool empty_col2) {
    out << indent << ind_newline << store_remember_font() << "</I></NOBR>" 
	<< ind_newline << "</TD></TR>";
    if ( ! empty_col2)
        out << "<TR><TD WIDTH=" << table_2c_first_col 
	    << "% NOWRAP></TD><TD ALIGN=LEFT VALIGN=TOP WIDTH="
	    << table_2c_second_col << "%>";
    out << outdent << ind_newline;
}

void two_cols_html_new_closing( ostream& out) {
    out << indent << ind_newline << store_remember_font() << "</I></NOBR>" 
	<< ind_newline << "</TD></TR>" << outdent << outdent;
}

void two_cols_html_end( ostream& out, bool empty_col2) {
    out << indent << ind_newline;
    if ( ! empty_col2)    
        out << "</TD></TR>" << ind_newline;
    out	<< "</TABLE><!2>" << outdent << outdent << ind_newline;
}


/* template declarations for functions and others */
/* ============================================== */
/* It depends on the layout (2 or 3 columns)      */

const char* handle_template_layout( ostream& out, 
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
	    if ( ! macroIsTrue( "\\ccTagRmTemplate")) {
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
		// Make the postprocessing glue the template keyword and list
		// together with the declaration.
		out << "<!GLUE>";
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
    while ( s != s_end && ( ! isalnum( *s)) && ( *s != '_') && ( *s != '|'))
        --s;
    // parse first identifier for the function_name
    while ( s != s_end && ( isalnum( *s) || *s == '_' || *s == '|'))
        --s;
    // check the possibilty that the function is a cast operator
    // like `operator int', or an idfier as operator like `operator new'
    // so check whether the next identifier is equal to `operator'
    p = s;
    while ( s != s_end && *s <= ' ')
        --s;
    while ( s != s_end && ( isalnum( *s) || *s == '_' || *s == '|'))
        --s;
    // check for conversion operators with template classes like
    // operator A<int>(...)
    const char* pp = signature;
    while ( isspace( *pp))
	pp++;
    if ( strncmp( pp, "operator", 8) == 0)
	s = pp - 1;
    else {
	if ( strncmp( s + 1, "operator", 8) != 0)
	    // it's not the cast operator, restore old position
	    s = p;
    }
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
	while ( s != s_end && ( isalnum( *s) || *s == '_' || *s == '|'))
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
    while ( s != s_end && ( ! isalnum( *s)) && ( *s != '_') && ( *s != '|'))
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
    while ( s != s_end && ( isalnum( *s) || *s == '_' || *s == '|'))
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
	while ( s != s_end && ( isalnum( *s) || *s == '_' || *s == '|'))
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
    if ( ! macroIsTrue( "\\ccTagRmConstRefPair"))
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

void remove_own_classname( char* s, const string& classname) {
    if ( ! macroIsTrue( "\\ccTagRmEigenClassName"))
        return;
    if ( classname.empty())
        return;
    // An easy algorithm (quadratic runtime worst case).
    // The classname has to match exactly (and is removed at most once).
    // If the result is an empty string, the classname is not removed.
    int l  = strlen( s) - classname.size();
    int l2 = classname.size();
    int i;
    bool is_empty = true;
    for ( i = 0; i <= l; i++) {
        if ( strncmp( s, classname.c_str(), l2) == 0)
	    break;
	if ( ! isspace( *s))
	    is_empty = false;
	s++;
    }
    if ( i <= l) {
        // Classname found. Test if it is not part of a larger idfier.
	if ( i < l && ( isalnum( s[l2]) || s[l2] == '_'))
	    return;
	if ( i > 0  && ( isalnum( s[-1]) || s[-1] == '_'))
	    return;
        // Classname found. Test non-empty result.
	char* p = s + l2;
	while ( *p) {
	    if ( ! isspace( *p))
		is_empty = false;
	    p++;
	}
	if ( ! is_empty) {
	    // Remove it.
	    for( ; l2 > 0; l2--)
		*s++ = ' ';
	}
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
    current_font = it_font;
    // Assume that the syntax is not malformed (error messages are
    // provided elsewhere, namely compiler or cc_manual.sty).
    // Non matching operators are printed in functional notation.
    if ( ! macroIsTrue( "\\ccTagOperatorLayout"))
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
	} else if (op[0] == '(' && op[1] == ')') {
	    infix   = "(";
	    postfix = ")";
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
	    infix     = template_class_name.c_str();
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


// Three column layout functions
// =====================================================

void format_function( bool method, const char* signature, 
		      bool is_empty_comment) {
    current_font = it_font;
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
        normal_operator     = (return_value != 0);
	conversion_operator = (return_value == 0);
    }

    double exp_size_ret = 0.0;
    double exp_size = 0.0;
    if ( return_value) {
        remove_const_ref_pair( return_value);
	exp_size_ret += estimate_html_size( return_value);
    }
    if ( conversion_operator)
	exp_size_ret += estimate_html_size( op_symbols);
 
    three_cols_html_begin( *current_ostream, 
			   exp_size_ret * stretch_factor > 
                           table_width*table_first_col/100.0);
    // ---------
    // index
    if ( !method && macroIsTrue( "\\lciIfHtmlIndex") 
	 && macroIsTrue( "\\lciIfHtmlRefIndex")) {
	char* formatted_return   = convert_fontified_ascii_to_html( 
	    return_value);
	char* formatted_scope    = convert_fontified_ascii_to_html( scope);
	char* formatted_function = convert_fontified_ascii_to_html( 
	    function_name);
	char* formatted_params   = convert_fontified_ascii_to_html(
	    parameter_list);
	char* formatted_rest     = convert_fontified_ascii_to_html( rest);
	// Make only the function name a hyperlink. Types could be substituted
	// according the rules.
	*index_stream << sort_key_function << '1';
	filter_for_index_comment( *index_stream, function_name);
	*index_stream << ' ';
	filter_for_index_comment( *index_stream, signature);
	*index_stream << "!><UL><LI><I>" 
		      << formatted_return << ' ' << formatted_scope
		      << "<A HREF=\"" << current_filename << "#Function_";
	filter_for_index_anchor( *index_stream, signature);
	*index_stream << "\">" << formatted_function << "</A>"
		      << "( " << formatted_params << ")" << formatted_rest
		      << "</I></UL>" << endl;
	*current_ostream << "<A NAME=\"Function_";
	filter_for_index_anchor( *current_ostream, signature);
	*current_ostream << "\"></A>" << endl;
	delete[] formatted_return;
	delete[] formatted_scope;
	delete[] formatted_function;
	delete[] formatted_params;
	delete[] formatted_rest;
    }
    // end index
    // ----------
   
    if ( return_value)
        print_ascii_to_html_spc( *current_ostream, return_value);
    if ( conversion_operator)
        print_ascii_to_html_spc( *current_ostream, op_symbols);

    // handle function body or operation signature
    // first, estimate size
    if ( method)
        exp_size = estimate_html_size( macroX( "\\ccPureVar")) + 
	           width_per_character;
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
    if (! macroIsTrue( "\\ccTagRmTrailingConst"))
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
		    *current_ostream,
		    exp_size_ret * stretch_factor > 
		        table_width * table_first_col / 100.0,
		    is_empty_comment || (exp_size * stretch_factor > 
                        table_width * table_second_col / 100.0)
		);
	    print_ascii_to_html_spc( *current_ostream, praefix);
	    *current_ostream << " ";
	    char* p = parameter_list;
	    if ( ! ignore_params) {
	        if ( method)
		    print_ascii_to_html_spc(*current_ostream, 
					    macroX( "\\ccPureVar"));
		else if (n) {
		    --n;
		    print_ascii_to_html_spc(*current_ostream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		}
	    }
	    *current_ostream << " ";
	    print_ascii_to_html_spc( *current_ostream, infix);
	    *current_ostream << " ";
	    if ( ! ignore_params && n) {
	        while ( n--) {
		    print_ascii_to_html_spc(*current_ostream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		    if ( n)
		        *current_ostream << ", ";
		}
	    }
	    print_ascii_to_html_spc( *current_ostream, postfix);
	}
    }
    if ( ! normal_operator || failed) {
	if ( scope)
	    exp_size += estimate_html_size( scope);
	exp_size += estimate_html_size( function_name);
	exp_size += 3 * width_per_character;  // for parameter list parantheses

	// then, do the printing
	three_cols_html_second(
		*current_ostream,
		exp_size_ret  * stretch_factor > 
		    table_width * table_first_col / 100.0,
		is_empty_comment || (exp_size * stretch_factor > 
                    table_width * table_second_col / 100.0)
	    );
	if ( conversion_operator) {
	    print_ascii_to_html_spc( *current_ostream, op_symbols);
	    *current_ostream << " ( ";
	    if ( ! definedMacro( "\\ccPureVar")) {
		printErrorMessage( VariableUsedError);
		*current_ostream << "*this";
	    } else
		print_ascii_to_html_spc( *current_ostream,
					 macroX( "\\ccPureVar"));
	    *current_ostream << ")";
	} else {
	    double dd_width = table_width * ( 1.0 - table_first_col / 100.0);
	    dd_width /= stretch_factor;
	    if ( exp_size > dd_width && parameter_list) {
		*current_ostream << store_remember_font();
	        *current_ostream << "<TABLE BORDER=0 CELLSPACING=0 "
		  "CELLPADDING=0><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP";
		if ( macroIsTrue( "\\ccLongParamLayout") ||
		     ! macroIsTrue( "\\ccAlternateThreeColumn"))
		    *current_ostream << " COLSPAN=2";
		*current_ostream << "><I>" << get_remember_font()
				 << ind_newline;
	    }
	    if ( method) {
		if ( ! definedMacro( "\\ccPureVar")) {
		    printErrorMessage( VariableUsedError);
		    *current_ostream << "*this";
		} else
		    print_ascii_to_html_spc(*current_ostream, 
					    macroX( "\\ccPureVar"));
		*current_ostream << '.';
	    }
	    if ( scope)
		print_ascii_to_html_spc( *current_ostream, scope);
	    print_ascii_to_html_spc( *current_ostream, function_name);
	    if ( parameter_list) {
		*current_ostream << " ( ";
		if ( exp_size > dd_width) {
		    *current_ostream << store_remember_font();
		    *current_ostream << "</I></TD>";
		    if ( macroIsTrue( "\\ccLongParamLayout") ||
			 ! macroIsTrue( "\\ccAlternateThreeColumn"))
			*current_ostream << "</TR><TR><TD WIDTH=" 
					<< table_long_param_indent 
					<< " NOWRAP></TD>";
		    *current_ostream << "<TD ALIGN=LEFT VALIGN=TOP "
		      "NOWRAP><I>";
		    *current_ostream << get_remember_font() << ind_newline;
		}
		char* p = parameter_list;
		while ( n--) {
		    print_ascii_to_html_spc( *current_ostream, p);
		    p += strlen( p) + 1;  // skip to next parameter
		    if ( n) {
			*current_ostream << ", ";
			if ( exp_size > dd_width)
			    *current_ostream << "<BR>" << ind_newline;
		    }
		}
		*current_ostream << ")";
		if ( exp_size > dd_width) {
		    *current_ostream << store_remember_font();
		    *current_ostream << "</I></TD></TR></TABLE>" <<ind_newline;
		}
	    } else
		*current_ostream << " ()";
	}
    }
    if (! macroIsTrue( "\\ccTagRmTrailingConst"))
        print_rest( *current_ostream, rest);
    three_cols_html_third( *current_ostream, 
			   exp_size  * stretch_factor> 
			        table_width * table_second_col / 100.0,
			   is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] function_name; 
    delete[] parameter_list; 
    three_cols_html_new_closing( *current_ostream);
}

void format_variable( const char* signature, 
		      bool is_empty_comment,
		      bool is_typedef = false) {
    current_font = it_font;
    char* return_value;
    char* scope;
    char* variable_name; 
    char* rest;    // possibly including assignment

    split_variable_declaration( signature, 
				return_value,
				scope,
				variable_name,
				rest);
    char* formatted_var = convert_fontified_ascii_to_html( variable_name);
    char* formatted_class = convert_fontified_ascii_to_html( template_class_name);

    double exp_size_ret = 0.0;
    double exp_size     = 0.0;
    if ( return_value) {
        remove_const_ref_pair( return_value);
	exp_size_ret += estimate_html_size( return_value);
    }
 
    three_cols_html_begin( *current_ostream, 
			   exp_size_ret * stretch_factor > 
			       table_width*table_first_col/100.0);
   
    if ( class_name.empty()) {
	if ( macroIsTrue( "\\lciIfHtmlLinks") && 
	     macroIsTrue( "\\lciIfHtmlRefLinks") && 
	     strlen(variable_name) > 1) {
	    // generate a substitution rule for hyperlinking
	    *anchor_stream << "[a-zA-Z0-9_]\"" << formatted_var
			   << "\"    { ECHO; }" << endl;
	    *anchor_stream << '"' << formatted_var
			   << "\"/{noCCchar}    { wrap_anchor( \""
			   << current_filename
			   << (is_typedef ? "#Typedef_" :  "#Var_" );
	    filter_for_index_anchor( *anchor_stream, variable_name);
	    *anchor_stream << "\", yytext); }" << endl;
	}
	if ( macroIsTrue( "\\lciIfHtmlIndex") && 
	     macroIsTrue( "\\lciIfHtmlRefIndex")) {
	    // index
	    *index_stream << (is_typedef ? sort_key_typedef :
			      sort_key_variable) << '1'; 
	    filter_for_index_comment( *index_stream, variable_name);
	    *index_stream << "!><UL><LI><I>" << formatted_var
			  << "</I></UL>" << endl;
	}
    } else  if ( macroIsTrue( "\\lciIfHtmlClassIndex") && 
		 macroIsTrue( "\\lciIfHtmlIndex")) {
	// index
	*index_stream << (is_typedef ? sort_key_typedef :
			  sort_key_variable) << '1'; 
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, variable_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename 
		      << (is_typedef ? "#Typedef_" :  "#Var_" );
	filter_for_index_anchor( *index_stream, variable_name);
	*index_stream << "\"><I>" << formatted_class << "::"
		      << variable_name << "</I></A></UL>" << endl;
    }
    if ( macroIsTrue( "\\lciIfHtmlLinks") || 
	 macroIsTrue( "\\lciIfHtmlIndex")) {
	*current_ostream << "<A NAME=\""
			<< (is_typedef ? "Typedef_" :  "Var_" );
	filter_for_index_anchor( *current_ostream, variable_name);
	*current_ostream << "\"></A>" << endl;
    }
    // end index

    if ( return_value)
        print_ascii_to_html_spc( *current_ostream, return_value);

    // handle function body or operation signature
    // first, estimate size
    if ( scope)
        exp_size += estimate_html_size( scope);
    exp_size += estimate_html_size( variable_name);
    if ( rest)
        exp_size += estimate_html_size( rest);

    // then, do the printing
    three_cols_html_second(
	    *current_ostream,
	    exp_size_ret * stretch_factor> table_width * table_first_col/100.0,
	    is_empty_comment || (exp_size  * stretch_factor> 
				 table_width * table_second_col / 100.0)
	);
    if ( scope)
        print_ascii_to_html_spc( *current_ostream, scope);
    print_ascii_to_html_spc( *current_ostream, variable_name);
    // *current_ostream << formatted_var;
    if ( rest) {
        *current_ostream << ' ';
        print_ascii_to_html_spc( *current_ostream, rest);
    }
    *current_ostream << ';';

    three_cols_html_third( *current_ostream, 
			   exp_size * stretch_factor >
			       table_width * table_second_col / 100.0,
			   is_empty_comment);
    delete[] return_value;
    delete[] scope;
    delete[] formatted_var; 
    delete[] variable_name; 
    delete[] rest; 
    three_cols_html_new_closing( *current_ostream);
}



// Two column layout functions
// =====================================================

void format_class_declaration( const char* signature) {
    current_font = it_font;
    char* return_value;
    char* scope;
    char* struct_name; 
    char* rest;

    split_variable_declaration( signature, 
				return_value,
				scope,
				struct_name,
				rest);

    char* formatted_struct = convert_fontified_ascii_to_html( struct_name);
    two_cols_html_begin( *current_ostream);

    if ( class_name.empty()) {
	if ( macroIsTrue( "\\lciIfHtmlLinks") && strlen(struct_name) > 1) {
	    // generate a substitution rule for hyperlinking
	    *anchor_stream << "[a-zA-Z0-9_]\"" << formatted_struct
			   << "\"    { ECHO; }" << endl;
	    *anchor_stream << '"' << formatted_struct
			   << "\"/{noCCchar}    { wrap_anchor( \""
			   << current_filename << "#Struct_";
	    filter_for_index_anchor( *anchor_stream, struct_name);
	    *anchor_stream << "\", yytext); }" << endl;
	}
	if ( macroIsTrue( "\\lciIfHtmlIndex")) {
	    // index
	    *index_stream << sort_key_struct << '1';
	    filter_for_index_comment( *index_stream, struct_name);
	    *index_stream << "!><UL><LI><I>" << formatted_struct
			  << "</I></UL>" << endl;
	}
    } else  if ( macroIsTrue( "\\lciIfHtmlClassIndex") && 
		 macroIsTrue( "\\lciIfHtmlIndex")) {
	// index
	*index_stream << sort_key_struct << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, struct_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Struct_";
	filter_for_index_anchor( *index_stream, struct_name);
	*index_stream << "\"><I>" << template_class_name << "::"
		      << struct_name << "</I></A></UL>" << endl;
    }
    if ( macroIsTrue( "\\lciIfHtmlLinks") || 
	 macroIsTrue( "\\lciIfHtmlIndex")) {
	*current_ostream << "<A NAME=\"Struct_";
	filter_for_index_anchor( *current_ostream, struct_name);
	*current_ostream << "\"></A>" << endl;
    }
    // end index

    print_ascii_to_html_spc( *current_ostream, signature);
    // *current_ostream << ';';

    delete[] return_value;
    delete[] scope;
    delete[] formatted_struct;
    delete[] struct_name;
    delete[] rest;

    two_cols_html_new_closing( *current_ostream);
}

void format_struct( const char* signature) {
    current_font = it_font;
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

    char* formatted_struct = convert_fontified_ascii_to_html( struct_name);
    two_cols_html_begin( *current_ostream);

    if ( class_name.empty()) {
	if ( macroIsTrue( "\\lciIfHtmlLinks") && strlen(struct_name) > 1) {
	    // generate a substitution rule for hyperlinking
	    *anchor_stream << "[a-zA-Z0-9_]\"" << formatted_struct
			   << "\"    { ECHO; }" << endl;
	    *anchor_stream << '"' << formatted_struct
			   << "\"/{noCCchar}    { wrap_anchor( \""
			   << current_filename << "#Struct_";
	    filter_for_index_anchor( *anchor_stream, struct_name);
	    *anchor_stream << "\", yytext); }" << endl;
	}
	if ( macroIsTrue( "\\lciIfHtmlIndex")) {
	    // index
	    *index_stream << sort_key_struct << '1';
	    filter_for_index_comment( *index_stream, struct_name);
	    *index_stream << "!><UL><LI><I>" << formatted_struct
			  << "</I></UL>" << endl;
	}
    } else  if ( macroIsTrue( "\\lciIfHtmlClassIndex") &&
		 macroIsTrue( "\\lciIfHtmlIndex")) {
	// index
	*index_stream << sort_key_struct << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, struct_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Struct_";
	filter_for_index_anchor( *index_stream, struct_name);
	*index_stream << "\"><I>" << template_class_name << "::"
		      << struct_name << "</I></A></UL>" << endl;
    }
    if ( macroIsTrue( "\\lciIfHtmlLinks") || 
	 macroIsTrue( "\\lciIfHtmlIndex")) {
	*current_ostream << "<A NAME=\"Struct_";
	filter_for_index_anchor( *current_ostream, struct_name);
	*current_ostream << "\"></A>" << endl;
    }
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
    if ( exp_size > table_width && parameter_list) {
      *current_ostream << store_remember_font();
      *current_ostream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	"<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP";
      //if ( tag_long_param_layout)
      //*current_ostream << " COLSPAN=2";
      *current_ostream << "><I>" << get_remember_font() << ind_newline;
    }
    if ( return_value) {
        print_ascii_to_html_spc( *current_ostream, return_value);
	*current_ostream << ' ';
    }
    if ( scope)
        print_ascii_to_html_spc( *current_ostream, scope);
    *current_ostream << formatted_struct;
    if ( parameter_list) {
        *current_ostream << " { ";
	if ( exp_size > table_width) {
	    *current_ostream << store_remember_font();
	    *current_ostream << "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
	    *current_ostream << get_remember_font() << ind_newline;
	}
	char* p = parameter_list;
	while ( n--) {
	    while ( *p && *p <= ' ')
	        ++p;
	    print_ascii_to_html_spc( *current_ostream, p);

	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_ostream << ", ";
		if ( exp_size > table_width)
		    *current_ostream << "<BR>" << ind_newline;
	    }
	}
        *current_ostream << "};";
	if ( exp_size > table_width) {
	    *current_ostream << store_remember_font();
	    *current_ostream << "</TD></TR></TABLE>" << ind_newline;
	}
    } else
        *current_ostream << ';';

    delete[] return_value;
    delete[] scope;
    delete[] formatted_struct;
    delete[] struct_name;
    delete[] parameter_list; 

    two_cols_html_new_closing( *current_ostream);
}

void format_nested_type( const char* nested_type_name) {
    current_font = it_font;
    char* formatted_type = convert_fontified_ascii_to_html( nested_type_name);
    char* formatted_class = convert_fontified_ascii_to_html( template_class_name);
    two_cols_html_begin( *current_ostream);

    if ( macroIsTrue( "\\lciIfHtmlClassIndex") && 
	 macroIsTrue( "\\lciIfHtmlIndex")) {
	// index
	*index_stream << sort_key_nested_type << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, nested_type_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Nested_type_" 
                      << nested_type_name
		      << "\"><I>";
	if ( ! template_class_name.empty())
	    *index_stream << formatted_class << "::";
	*index_stream << nested_type_name << "</I></A></UL>" << endl;
	*current_ostream << "<A NAME=\"Nested_type_" << nested_type_name 
			<< "\"></A>" << endl;
	// end index
    }

    // first, estimate size
    double exp_size =   estimate_html_size( template_class_name)
                      + estimate_html_size( nested_type_name) 
                      + 2 * width_per_character;
    exp_size *= stretch_factor;

    print_ascii_to_html_spc( *current_ostream, template_class_name);
    *current_ostream << "::";
    print_ascii_to_html_spc( *current_ostream, nested_type_name);

    delete[] formatted_type;

    two_cols_html_new_closing( *current_ostream);
}

// NOTE:  when an enum is formatted inside a ccRefEnum environment
// it currently always looks like it's inside a class (i.e. class_name
// is not empty), so the enum tags will NOT be linked or indexed 
// automatically.  This could be fixed by setting a flag of some sort
// to indicate that it's not really a class or perhaps finding a different
// way to make the classname and templateclassname known internally instead
// of using a \begin{lciClass} in \lciRefDeclContX (in cc_manual.sty)
void format_enum( const char* signature) {
    current_font = it_font;
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

    char* formatted_enum = convert_fontified_ascii_to_html( enum_name);

    two_cols_html_begin( *current_ostream);

    if ( class_name.empty()) {
	if ( macroIsTrue( "\\lciIfHtmlLinks")  && strlen( enum_name) > 1) {
	    // generate a substitution rule for hyperlinking
	    *anchor_stream << "[a-zA-Z0-9_]\"" << formatted_enum
			   << "\"    { ECHO; }" << endl;
	    *anchor_stream << '"' << formatted_enum
			   << "\"/{noCCchar}    { wrap_anchor( \""
			   << current_filename << "#Enum_";
	    filter_for_index_anchor( *anchor_stream, enum_name);
	    *anchor_stream << "\", yytext); }" << endl;
	}
	if ( macroIsTrue( "\\lciIfHtmlIndex")) {
	    // index
	    *index_stream << sort_key_enum << '1';
	    filter_for_index_comment( *index_stream, enum_name);
	    *index_stream << "!><UL><LI><I>" << formatted_enum
			  << "</I></UL>" << endl;
	}
    } else  if ( macroIsTrue( "\\lciIfHtmlClassIndex") &&
		 macroIsTrue( "\\lciIfHtmlIndex")) {
	// index
	*index_stream << sort_key_enum << '1';
	filter_for_index_comment( *index_stream, template_class_name);
	*index_stream << "::";
	filter_for_index_comment( *index_stream, enum_name);
	*index_stream << "!><UL><LI><A HREF=\""
		      << current_filename << "#Enum_";
	filter_for_index_anchor( *index_stream, enum_name);
	*index_stream << "\"><I>" << template_class_name << "::"
		      << enum_name << "</I></A></UL>" << endl;
    }
    if ( macroIsTrue( "\\lciIfHtmlLinks") || 
	 macroIsTrue( "\\lciIfHtmlIndex")) {
	*current_ostream << "<A NAME=\"Enum_";
	filter_for_index_anchor( *current_ostream, enum_name);
	*current_ostream << "\"></A>" << endl;
    }
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
    if ( exp_size > table_width && parameter_list) {
      *current_ostream << store_remember_font();
      *current_ostream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	"<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP";
      //if ( tag_long_param_layout)
      //*current_ostream << " COLSPAN=2";
      *current_ostream << "><I>" << get_remember_font() << ind_newline;
    }
    if ( return_value) {
        print_ascii_to_html_spc( *current_ostream, return_value);
	*current_ostream << ' ';
    }
    if ( scope)
        print_ascii_to_html_spc( *current_ostream, scope);
    *current_ostream << formatted_enum;
    if ( parameter_list) {
        *current_ostream << " { ";
	if ( exp_size > table_width) {
	    *current_ostream << store_remember_font();
	    *current_ostream << "</TD><TD ALIGN=LEFT VALIGN=TOP NOWRAP>";
	    *current_ostream << get_remember_font() << ind_newline;
	}
	char* p = parameter_list;
	while ( n--) {
	    while ( *p && *p <= ' ')
	        ++p;
	    print_ascii_to_html_spc( *current_ostream, p);

	    if ( class_name.empty()) {
		if ( macroIsTrue( "\\lciIfHtmlIndex")) {
		    // index: print enum tags with (possible) initializers
		    *index_stream << sort_key_enum_tags  << '1' << p 
				  << "!><UL><LI><I>";
		    print_ascii_to_html_spc( *index_stream, p);
		    *index_stream << "</I></UL>" << endl;
		}
		if ( macroIsTrue( "\\lciIfHtmlLinks")) {
		    // generate a substitution rule for hyperlinking
		    // Here, the initializer has to be suppressed
		    char* q = p;
		    while( *q && *q != '=')
			++q;
		    while ( q>p && q[-1] <= ' ')
			--q;
		    char c_tmp = *q;
		    *q = '\0';
		    if ( strlen( p) > 1) {
			char *tmp_param = convert_fontified_ascii_to_html( p);
			*anchor_stream << "[a-zA-Z0-9_]\"" << tmp_param
				  << "\"    { ECHO; }" << endl;
			*anchor_stream << '"' << tmp_param
				  << "\"/{noCCchar}    { wrap_anchor( \""
				  << current_filename << "#Enum_";
			filter_for_index_anchor( *anchor_stream, enum_name);
			*anchor_stream << "\", yytext); }" << endl;
			delete[] tmp_param;
		    }
		    *q = c_tmp; // restore initializer
		}
            } else  if ( macroIsTrue( "\\lciIfHtmlClassIndex") && 
			 macroIsTrue( "\\lciIfHtmlIndex")) {
		// index: print enum tags with their (possible) initializers
		*index_stream << sort_key_enum_tags  << '1' << p 
			      << "!><UL><LI><I>";
		print_ascii_to_html_spc( *index_stream, p);
		*index_stream << "</I></UL>" << endl;

		*index_stream << sort_key_enum_tags << '1';
		filter_for_index_comment( *index_stream, template_class_name);
		*index_stream << "::";
		filter_for_index_comment( *index_stream, p);
		*index_stream << "!><UL><LI><A HREF=\""
			      << current_filename << "#Enum_";
		filter_for_index_anchor( *index_stream, enum_name);
		*index_stream << "\"><I>" << template_class_name << "::"
			      << p << "</I></A></UL>" << endl;
	    }

	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_ostream << ", ";
		if ( exp_size > table_width)
		    *current_ostream << "<BR>" << ind_newline;
	    }
	}
        *current_ostream << "};";
	if ( exp_size > table_width) {
	    *current_ostream << store_remember_font();
	    *current_ostream << "</TD></TR></TABLE>" << ind_newline;
	}
    } else
        *current_ostream << ';';

    delete[] return_value;
    delete[] scope;
    delete[] formatted_enum;
    delete[] enum_name;
    delete[] parameter_list; 

    two_cols_html_new_closing( *current_ostream);
}

void format_constructor( const char* signature) {
    current_font = it_font;
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

    two_cols_html_begin( *current_ostream);
    // first, estimate size
    double exp_size = 1.0;
    // Of no use here!
    // if ( scope)
    //	exp_size += estimate_html_size( scope);
    exp_size += estimate_html_size( macroX( "\\ccPureVar")) 
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
    if ( exp_size > table_width && parameter_list) {
      *current_ostream << store_remember_font();
      *current_ostream << "<TABLE BORDER=0 CELLSPACING=0 CELLPADDING=0>"
	"<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP";
      if ( macroIsTrue( "\\ccLongParamLayout") ||
	   ! macroIsTrue( "\\ccAlternateThreeColumn"))
	  *current_ostream << " COLSPAN=2";
      *current_ostream << "><I>" << get_remember_font() << ind_newline;
    }
    // Of no use here!
    // if ( scope)
    //    print_ascii_to_html_spc( *current_ostream, scope);
    print_ascii_to_html_spc( *current_ostream, template_class_name);
    *current_ostream << ' ';  // More space would be a good idea here.
    print_ascii_to_html_spc( *current_ostream, macroX( "\\ccPureVar"));

    if ( parameter_list) {
        *current_ostream << " ( ";
	if ( exp_size > table_width) {
	  *current_ostream << store_remember_font();
	  *current_ostream << "</I></TD>";
	  if ( macroIsTrue( "\\ccLongParamLayout") ||
	       ! macroIsTrue( "\\ccAlternateThreeColumn"))
	      *current_ostream << "</TR><TR><TD WIDTH=" 
			      << table_long_param_indent 
			      << " NOWRAP></TD>";
	  *current_ostream << "<TD ALIGN=LEFT VALIGN=TOP "
	    "NOWRAP><I>";
	  *current_ostream << get_remember_font() << ind_newline;
	}
	char* p = parameter_list;
	while ( n--) {
	    print_ascii_to_html_spc( *current_ostream, p);
	    p += strlen( p) + 1;  // skip to next parameter
	    if ( n) {
	        *current_ostream << ", ";
		if ( exp_size > table_width)
		    *current_ostream << "<BR>" << ind_newline;
	    }
	}
        *current_ostream << ");";
	if ( exp_size > table_width) {
	    *current_ostream << store_remember_font();
	    *current_ostream << "</I></TD></TR></TABLE>" << ind_newline;
	}
    } else
        *current_ostream << ';';

    delete[] return_value;
    delete[] scope;
    delete[] function_name; 
    delete[] parameter_list; 

    two_cols_html_new_closing( *current_ostream);
}

// Toplevel handleFunctions
// =====================================================

void handle_two_column_layout( char key, const char* decl) {
    if ( current_ostream) {
        decl = handle_template_layout( *current_ostream, decl, false);
	switch ( key) {
	case 'A':
	    format_class_declaration( decl);
	    break;
	case 'B':
	    format_struct( decl);
	    break;
	case 'C':
	    format_nested_type( decl);
	    break;
	case 'D':
	    format_enum( decl);
	    break;
	case 'E':
	    format_constructor( decl);
	    break;
	default:
	    printErrorMessage( UnknownKeyError);
	}
    }
}

void handle_three_column_layout( char key, const char* decl, bool empty) {
    if ( current_ostream) {
        decl = handle_template_layout( *current_ostream, decl, true);
	switch ( key) {
	case 'L':  // member function
	    format_function( true, decl, empty);
	    break;
	case 'M':  // function
	    format_function( false, decl, empty);
	    break;
	case 'N':
	    format_variable( decl, empty);
	    break;
	case 'O': // typedef
	    format_variable( decl, empty, true);
	    break;
	default:
	    printErrorMessage( UnknownKeyError);
	}
    }
}


// EOF //

