#ifndef REGEX_DICTIONARY_H
#define REGEX_DICTIONARY_H

#include <mstring.h>

void          regex_register( const string& name, const string &regex );
bool          regex_does_match( const string& name, const string& text );
const string& regex_get_submatch( unsigned int num );

#endif
