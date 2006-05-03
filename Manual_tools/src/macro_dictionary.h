/**************************************************************************

  macro_dictionary.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Dictionary of TeX macro definitions.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$

**************************************************************************/

#ifndef MACRO_DICTIONARY_H
#define MACRO_DICTIONARY_H 1

#include <mstring.h>

bool is_active_char( char c);

typedef string (*ExpandFunction)( const string& macro,
				  string parameters[],
				  size_t n_parameters,
				  size_t n_options);

struct Macro_item {
    string  filename;
    size_t  line;
    string  body;
    size_t  n_param;
    size_t  n_opt_at_end;
    ExpandFunction fct;

    Macro_item() : n_param(0), n_opt_at_end(0) {}
    Macro_item(	const string&  file,
		size_t         ln,
		const string&  bdy,
		size_t         n_par = 0)
	: filename(file), line(ln), body(bdy), n_param(n_par),
	  n_opt_at_end(0), fct(0)
    {}
    Macro_item(	const string&  file,
		size_t         ln,
		ExpandFunction f,
		size_t         n_par = 0)
	: filename(file), line(ln), n_param(n_par), n_opt_at_end(0), fct(f)
    {}
};

void          pushMacroScope();
void          popMacroScope();

void          insertMacro( const string& macro,
			   const string& filename,
			   size_t        line,
			   const string& body,
			   size_t        n_param = 0);

void          insertInternalMacro( const string&  macro,
				   ExpandFunction fct,
				   size_t         n_param = 0);

void          insertInternalMacro( const string&  macro,
				   const string&  body,
				   size_t         n_param = 0);

void          insertGlobalMacro( const string& macro,
				 const string& filename,
				 size_t        line,
				 const string& body,
				 size_t        n_param = 0);

void          insertInternalGlobalMacro( const string&  macro,
					 ExpandFunction fct,
					 size_t         n_param = 0);

void          insertInternalGlobalMacro( const string&  macro,
					 const string&  body,
					 size_t         n_param = 0);

const Macro_item& fetchMacro(  const string& macro);

inline
const string& fetchMacroBody(  const string& macro) {
    return fetchMacro( macro).body;
}

inline
const string& macroX(  const string& macro) {
    return fetchMacroBody( macro);
}

void eraseMacro( const string& macro);

string expandFirstMacro( string body, bool expand_only_once = false);

string expandMacro( const string& macro,
		    const Macro_item& item,
		    string parameters[],
		    size_t n_parameters,
		    size_t n_options);


bool definedMacro(  const string& macro);

bool macroIsTrue( const string& macro);

void checkMacroOptEnd(  const string& macro);

#endif // MACRO_DICTIONARY_H //
// EOF //

