#include <iostream>
#include <mstring.h>
#include <cassert>
#include "regex++.h"

using namespace std;

const unsigned int maxgroups = 20;

vector<string>
Regular_expression::submatch(maxgroups);


Regular_expression::Regular_expression( const char* rx )
{ init( rx ); }

Regular_expression::Regular_expression( const string& rx )
{ init( rx.c_str() ); }

void
Regular_expression::init( const char* rx ) {
  assert( rx != NULL );
  if( regcomp( &regex, rx, REG_EXTENDED ) != 0 ) {
    std::cout << "!! Error: regex compilation error" << std::endl;
    exit(1);
  }
	regmatch = ( regmatch_t* ) malloc ( maxgroups * sizeof ( regmatch_t ) );
  assert( regmatch != NULL );
}

Regular_expression::~Regular_expression() {
  if( regmatch != NULL ) {
    regfree ( &regex );
    free ( regmatch );
  }
}

bool
Regular_expression::match( const string& s ) {
  assert( regmatch != NULL );
	if( regexec( &regex, s.c_str(), maxgroups, regmatch, 0 ) == 0 )	{
		int so, eo;
		for ( unsigned int i = 0; i < maxgroups; i++ ) {
			so = regmatch[ i ].rm_so;
			eo = regmatch[ i ].rm_eo;
			if ( so != -1 )
				submatch[ i ] = s.substr ( so, eo - so );
			else
				submatch[ i ] = "";
		}
		return true;
	}
  return false;
}

bool
Regular_expression::match( const char *s )
{ return Regular_expression::match( string( s ) ); }

const string&
Regular_expression::get_match( unsigned int num ) {
  static string emptystring;
  if( num >= maxgroups ) {
    std::cerr << "!! warning: accessing bad regex submatch #" << num << std::endl;
    return emptystring;
  }

  return submatch[num];
}



