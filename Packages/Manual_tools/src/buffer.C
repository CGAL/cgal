/**************************************************************************
 
  database.C
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Datastructures for the database extracted from
              the TeX and C++ code mixed files.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

// Other modules needed:
// =====================
#define ADT_INCLUDE_CC 1
#include <database.h>
#include <ADT/include_cc.h>


// Substitute old style malloc, realloc, strdup ...
// ================================================
char* newstr( const char* src) {
    ADT_Assert( src);
    if ( ! src)
        return 0;
    char* s = new char[ strlen( src) + 1];
    strcpy( s, src);
    return s;
}


// Class definitions:
// ==============================================
int printComment( ostream &out, const Text& T, bool leadingLine, bool HTML) {
    InListFIter< TextToken> words( (Text&)T);
    int  width = MaxTextWidth - indentation_number() - 3;
    int  w     = width;
    int  state = 0;     // 0 = start, 1 = after token, 2 = spaces
                        // 3 = one newline (and spaces), 4 = newlines
    ForAll( words) {
	switch ( state) {
	case 0:
	    if ( !words->isSpace) {
		if ( leadingLine)
		    out << endl;
		out << indNewline << (HTML ? "" : "// ") << words->string;
		w -= words->len;
		state = 1;
	    }
	    break;
	case 1:
	    if ( !words->isSpace || words->len > 0) {
		if ( !words->isSpace) {
		    if ((words->len > w) && (w != width)) {
			w = width;
			out << indNewline << (HTML ? "" : "// ");
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
			out << indNewline << (HTML ? "" : "// ");
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
			out << indNewline << (HTML ? "" : "// ");
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
		    if (HTML)
		        out << indNewline << "<P>";
		    out << indNewline << (HTML ? "" : "// ");
		    out << indNewline << (HTML ? "" : "// ");
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

int printTrueComment( ostream &out, const Text& T, bool leadingLine) {
    InListFIter< TextToken> words( (Text&)T);
    int state  = 0;  // 0 = start, 1 = spaces at start, 2 = at least one token
    int spaces = 0;
    ForAll( words) {
	switch ( state) {
	case 0:
	    if ( !words->isSpace) {
		if ( leadingLine)
		    out << endl;
		out << indNewline << "// " << words->string;
		state = 1;
	    } else {
		if ( words->string[0] != '\n') {
		    state  = 1;
		    spaces = 1;
		}
	    }
	    break;
	case 1:
	    if ( !words->isSpace) {
		if ( leadingLine)
		    out << endl;
		out << indNewline << "// ";
		while ( spaces--)
		    out << ' ';
		out << words->string;
	    } else {
		if ( words->string[0] == '\n')
		    state = 0;
		else
		    spaces++;
	    }
	    break;
	case 2:
	    if ( words->isSpace && ( words->string[0] == '\n'))
		out << indNewline << "// "; 
	    else
		out << words->string;
	    break;
	}
    }
    return state;
}


ostream&
operator<< (ostream& out, const TextToken& t) {
    out << t.isSpace << ", ";
    printString( out, t.string, t.len);
    return out;
}

istream&
operator>> (istream& in, TextToken& t) {
    char c;
    int  e;
    in >> e >> c;  // c == ','
    t.isSpace = bool( e);
    delete[] t.string;
    t.len = scanString( in, t.string);
    return in;
}

ostream&
operator<< (ostream& out, const Declaration& d) {
    out << "Declaration[" << indent << indNewline;
    out << d.type << indNewline;
    printString( out, d.returnType);
    out << indNewline;
    printString( out, d.name);
    out << indNewline;
    printString( out, d.templateParams);
    out << indNewline;
    printString( out, d.parameters);
    out << indNewline;
    out << d.comment << indNewline;
    out << d.spec << outdent << indNewline;
    out << ']';
    return out;
}

istream&
operator>> (istream& in, Declaration& ) {
    return in;
}

ostream&
operator<< (ostream& out, const Specification& s) {
    out << "Specification[" << indent << indNewline;
    out << s.type    << indNewline;
    out << s.decl    << indNewline;
    out << s.comment << outdent << indNewline;
    out << ']';
    return out;
}

istream&
operator>> (istream& in, Specification& ) {
    return in;
}

// Auxiliary function definitions
// ==============================================
/* ...
    char* catString( const char* r, const char* s, int l, int m) {
    char *p;
    if ( l < 0)
	l = strlen( r);
    if ( m < 0)
	m = strlen( s);
    p = (char *)malloc( l + m + 1);
    ADT_Assert( p != NULL);
    memcpy( p, r, l);
    memcpy( p+l, s, m);
    p[ l+m] = 0;
    return( p);
    }
... */

void printString( ostream &out, const char *s, int len) {
    if ( len < 0)
	len = strlen( s);
    out << len;
    if ( len > 0)
        out << ", " << s;
}

int scanString(  istream &in,  char*& s) {
    int len;
    int i;
    in >> len;
    if ( len > 0) {
	s = new char[ len + 1];
	ADT_Assert( s != NULL);
	in.get( s[0]);  // the comma
	in.get( s[0]);  // the space
	for ( i=0; i<len; i++)
	    in.get( s[i]);
	s[i] = 0;
    } else {
	s = NULL;
    }
    return len;
}



// EOF //

