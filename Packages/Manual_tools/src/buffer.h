/**************************************************************************
 
  database.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Datastructures for the database extracted from
              the TeX and C++ code mixed files.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_DATABASE)
#define MODULE_DATABASE 1

#define MaxTextWidth 74

// Other modules needed:
// ==============================================
#include <stdlib.h>
#include <ADT/adtbase.h>
#undef LK_RestrictedOverloading
#define LK_RestrictedOverloading 0
#include <ADT/list.h>
#include <string.h>
#include <ctype.h>

// Class declarations:
// ==============================================
// Top Level Structure:
// typedef InList<Specification> SpecificationList;  // typedef is postponed.
struct Specification;

// Parts of a specification:
struct Declaration;
class  TextToken;
// typedef InList< TextToken> Text;  // typedef is postponed.


// Substitute old style malloc, realloc, strdup ...
// ================================================

char* renew( char* old, size_t old_size, size_t new_size);
char* newstr( const char* src);


// Class specifications:
// ==============================================

// Auxiliary class:
// ------------------------------------
class Buffer {
    size_t fibo1;       // fibonacci numbers to increase buffer in
    size_t fibo2;       // the case of an overflow.
    size_t len;
    char*  data;
public:
    Buffer()  {
	fibo1 = 610;
	fibo2 = 987;
	len   = 0;
	data  = new char[ fibo2];
	*data = '\0';
	// Assert( data != NULL);
    }
    ~Buffer() {
        ADT_Assert( fibo2 > fibo1);  // SGI complains in the test_suit
        ADT_Assert( fibo2 > len);   // Gnu complains in the test_suit
        ADT_Assert( data[ len] == '\0');
        ADT_Assert( strlen( data) == len);
	delete[] data;
    }
    Buffer( const Buffer& t) {
        ADT_Assert( t.data[ t.len] == '\0');
	fibo1 = t.fibo1;
	fibo2 = t.fibo2;
	len   = t.len;
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
	data  = newstr( t.data); // ==============  ?? FALSCH ERROR !!!!
	cerr << "Error: Buffer copy constructor called." << endl;
	exit( 6);
    }
    Buffer& operator=( const Buffer& t) {
	if( this != &t) { // beware of this=t; see Stroustrup pp.238
	    delete[] data;
            t.data[ t.len] = '\0';
	    fibo1 = t.fibo1;
	    fibo2 = t.fibo2;
	    len   = t.len;
	    data  = newstr( t.data);// ==============  ?? FALSCH ERROR !!!!
	} else {
	    fibo1 = 610;
	    fibo2 = 987;
	    len   = 0;
	    data  = new char[ fibo2];
	}
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
	cerr << "Error: Buffer assignment operator called." << endl;
	exit( 6);
	return *this;
    }
    void add( const char* s, int l = -1) {
        ADT_Assert( data[ len] == '\0');
	int tmp;
	if ( l < 0)
	    l = strlen( s);
	// ADT_Assert( l == (int)strlen( s));  // this is definitly wrong!
	if ( len + l >= fibo2) { // increase buffer
	    do {
		tmp    = fibo2;
		fibo2 += fibo1;
		fibo1  = tmp;
	    } while ( len + l >= fibo2);  // increase buffer
	    data = renew( data, fibo1, fibo2);
	    // Assert( data != NULL);
	}
	memcpy( data + len, s, l);
	len += l;
	data[ len] = '\0';
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
    }
    void add( const Buffer& buf) {
        add( buf.string(), buf.length());
    }
    void add( const Buffer* buf) {
        add( buf->string(), buf->length());
    }
    void add( char c) {
        ADT_Assert( data[ len] == '\0');
	if ( len + 1 >= fibo2) { // increase buffer
	    int tmp = fibo2;
	    fibo2  += fibo1;
	    fibo1   = tmp;
	    data = renew( data, fibo1, fibo2);
	    //	    Assert( data != NULL);
	}
	data[ len++] = c;
	data[ len]   = '\0';
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
    }
    void prepend( const char* s, int l = -1) {
	ADT_Assert( s);
        ADT_Assert( data[ len] == '\0');
	int tmp;
	if ( l < 0)
	    l = strlen( s);
	ADT_Assert( l == (int)strlen( s));
	if ( len + l >= fibo2) { // increase buffer
	    do {
		tmp    = fibo2;
		fibo2 += fibo1;
		fibo1  = tmp;
	    } while ( len + l >= fibo2);  // increase buffer
	    data = renew( data, fibo1, fibo2);
	    // Assert( data != NULL);
	}
	// relocate current text
	int i;
	for ( i = len-1; i >= 0; i--)
	    data[i+l] = data[ i];
	memcpy( data, s, l);
	len += l;
	data[ len]   = '\0';
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
    }
    void prepend( char c) {
        ADT_Assert( data[ len] == '\0');
	if ( len + 1 >= fibo2) { // increase buffer
	    int tmp = fibo2;
	    fibo2  += fibo1;
	    fibo1   = tmp;
	    data = renew( data, fibo1, fibo2);
	    //	    Assert( data != NULL);
	}
	// relocate current text
	int i;
	for ( i = len; i > 0; i--)
	    data[i] = data[ i-1];
	data[ 0] = c;
	len++;
	data[ len]   = '\0';
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
    }
    void cut( int i) {
        ADT_Assert( data[ len] == '\0');
        ADT_Assert( i >= 0);
        if ( (int)len >= i) {
	    len -= i;
	}
	data[ len]   = '\0';
        ADT_Assert( fibo2 > fibo1);
        ADT_Assert( fibo2 > len);
    }
    const char* string() const {
        ADT_Assert( data[ len] == '\0');
	ADT_Assert( len < fibo2);
	return data;
    }
    void set( int i, char c) {
	ADT_Assert( c != 0);
	ADT_Assert( i >= 0);
        ADT_Assert( size_t(i) < len);
	data[i] = c;
    }
    int length() const {
        ADT_Assert( data[ len] == '\0');
	ADT_Assert( len < fibo2);
	return (int)len;
    }
    void flush() {
        ADT_Assert( data[ len] == '\0');
	ADT_Assert( len < fibo2);
	len = 0;
	*data = '\0';
    }
    void capitalize() {
	for ( int i = 0; i < len; i++) 
	    data[i] = toupper( data[i]);
    }
};

class TextToken : public ListLink{
public:
    bool  isSpace;
    int   len;
    char* string;
    TextToken() {
	isSpace = false;
	len     = 0;
	string  = 0;
    }
    TextToken( const char* s, int l = -1, bool space = false) {
        isSpace = space;
	if ( s) {
	    if ( l < 0)
		len = strlen( s);
	    else
		len = l;
	    string  = newstr( s);
	} else {
	    string  = newstr( "");
	    len     = 0;
	}
    }
    TextToken( const TextToken& t) {
        isSpace = t.isSpace;
	len     = t.len;
	string  = newstr( t.string);
    }
    ~TextToken() {
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	delete[] string;
    }
    TextToken& operator=( const TextToken& t) {
	if( this != &t) { // beware of this=t; see Stroustrup pp.238
	    delete[] string;
	    isSpace = t.isSpace;
	    len     = t.len;
	    string  = newstr( t.string);
	}
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	return *this;
    }
    TextToken& add( const TextToken& t) {
	if ( t.len > 0) {
	    string = renew( string, len + 1, len + t.len + 1);
	    memcpy( string + len, t.string, t.len);
	    len += t.len;
	    string[len] = '\0';
	    isSpace = isSpace && t.isSpace;
	}
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	return *this;
    }
    TextToken& add( char c) {
        string = renew( string, len + 1, len + 2);
	string[ len] = c;
	len ++;
	string[ len] = '\0';
	isSpace = isSpace && (( c == ' ') || ( c == '\t'));
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	return *this;
    }
    TextToken& prepend( const TextToken& t) {
	if ( t.len > 0) {
	    int i;
	    string = renew( string, len + 1, len + t.len + 1);
	    for( i=len; i>=0; i--)  // 0 will also be copied
		string[ i + t.len] = string[ i];
	    memcpy( string, t.string, t.len);
	    len += t.len;
	    isSpace = isSpace && t.isSpace;
	}
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	return *this;
    }
    TextToken& prepend( char c) {
	int i;
        string = renew( string, len + 1, len + 2);
	for( i=len; i>=0; i--)  // 0 will also be copied
	    string[ i + 1] = string[ i];
	*string = c;
	len ++;
	isSpace = isSpace && (( c == ' ') || ( c == '\t'));
        ADT_Assert( (! string && len == 0) || (int)strlen( string) == len);
        ADT_Assert( ! string || string[ len] == '\0');
	return *this;
    }
    friend ostream& operator<< (ostream& out, const TextToken& t);
    friend istream& operator>> (istream& in, TextToken& t);
};

//class   InList< TextToken>;

typedef InList< TextToken> Text;

int printComment(     ostream& out, const Text& T, 
                      bool leadingLine = false, bool HTML = false);
int printTrueComment( ostream& out, const Text& T, bool leadingLine = false);



void printString( ostream &out, const char*  s, int len = -1);
int  scanString(  istream &in,        char*& s);


#endif // MODULE_DATABASE //
