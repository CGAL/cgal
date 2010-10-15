/**************************************************************************
 
  basic.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Assertions, indented output, etc.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_BASIC)
#define MODULE_BASIC 1

#include <cstdlib>
#include <iostream>
#include <string.h>

using namespace std;

// Own assertion macro
// ================================================

void cc_assertion_error( char *cond, char *fname, int line);

#if !NDEBUG
#define CC_Assert( cond) {\
    if (!( cond)) \
        cc_assertion_error( #cond, __FILE__, __LINE__); \
    }
#else // NDEBUG
#define CC_Assert( cond) (void)0
#endif // NDEBUG


// Indented output, using global indentation counter
// ====================================================

ostream& indent(       ostream&);
ostream& outdent(      ostream&);
ostream& ind_newline(  ostream&);
ostream& ind_space(    ostream&);
int indentation_number();

// Substitute old style malloc, realloc, strdup ...
// ================================================

char* renew( char* old, size_t old_size, size_t new_size);
char* newstr( const char* src);

int execute_shell_command( const string& cmd, std::ostream& out, std::ostream& err );


#endif // MODULE_BASIC //
