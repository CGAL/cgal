/**************************************************************************
 
  basic.C
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Assertions, indented output, etc.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <basic.h>
#include <stdlib.h>
#include <string.h>

// Own assertion macro
// ================================================

void cc_assertion_error( char *cond, char *fname, int line) {
    cerr << endl << "fatal error: assertion '" << cond << "' failed in line " 
	 << line << " of file '" << fname << "'." << endl;
    abort();
}

// Indented output, using global indentation counter
// ====================================================

class Output_indentation {
    int n;
    int step;
public:
    Output_indentation( int stepSize = 4) : n(0), step( stepSize) {}
    void operator++(){ n += step;}
    void operator--(){ n -= step;}
    int indentation() { return n; }
};

ostream& operator<< (ostream& out, Output_indentation& ind) {
    const char* s = "                                "; // 32 spaces !!!
    int n = ind.indentation();
    while ( n > 32) { // 32 !!!
	n -= 32;
	out << s;
    }
    if ( n > 0)
	out << (s + 32 - n);
    return out;
}

static Output_indentation output_indentation;

ostream& indent( ostream& out){
    ++ output_indentation;
    return out;
}

ostream& outdent( ostream& out){
    -- output_indentation;
    return out;
}

ostream& ind_newline( ostream& out){
    out << endl << output_indentation;
    return out;
}

ostream& ind_space( ostream& out){
    out << output_indentation;
    return out;
}

int indentation_number() {
    return output_indentation.indentation();
}

// Substitute old style malloc, realloc, strdup ...
// ================================================
char* renew( char* old, size_t old_size, size_t new_size) {
    CC_Assert( old_size == 0 || old != 0);
    CC_Assert( new_size > old_size);
    char* cpy = new char[ new_size];
    if ( old && old_size > 0) {
	size_t min = ( old_size < new_size ? old_size : new_size);
	memcpy( cpy, old, min);
	delete[] old;
    }
    return cpy;
}

char* newstr( const char* src) {
    CC_Assert( src);
    if ( ! src)
        return 0;
    char* s = new char[ strlen( src) + 1];
    strcpy( s, src);
    return s;
}

// EOF //
