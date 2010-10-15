/**************************************************************************
 
  buffer.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : List of buffers (similar to strings) keeping the text
              parsed from the TeX input files.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_BUFFER)
#define MODULE_BUFFER 1


// Other modules needed:
// ==============================================
#include <stdlib.h>
#include <list>
#include <vector>
#include <basic.h>


// Class specifications:
// ==============================================
class Buffer {
    typedef std::vector<char>::const_iterator const_iterator;
    std::vector<char> buf;
public:
    Buffer() : buf( 1, '\0') {}
    Buffer( char c) : buf( 2, '\0') { *(buf.begin()) = c; }
    Buffer( const char* s, int l = -1) {
	if (s)
	    buf.insert( buf.begin(), s, s + 1 + ((l<0) ? strlen(s) : l));
	else
	    buf.push_back( '\0');
    }
    void add( const char* s, int l = -1) {
	if ( l < 0)
	    l = strlen( s);
	buf.insert( buf.end() - 1, s, s + l);
    }
    void add( const Buffer& buf) { add( buf.begin(), buf.size() - 1); }
    void add( char c) {
	buf.pop_back();
	buf.push_back(c);
	buf.push_back( '\0');
    }
    void prepend( const char* s, int l = -1) {
	if ( l == 0 || s == 0)
	    return;
	if ( l < 0)
	    l = strlen( s);
	buf.insert( buf.begin(), s, s + l);
    }
    void        prepend( const Buffer& t) { prepend( t.begin(), t.size() - 1);}
    void        prepend( char c)     { buf.insert( buf.begin(), c); }
    char*       begin()              { return &*(buf.begin()); }
    const char* begin()        const { return &*(buf.begin()); }
    char*       end()                { return &*(buf.end()); }
    const char* end()          const { return &*(buf.end()); }
    void        set( int i, char c)  { buf[i] = c; }
    size_t      size()         const { return buf.size(); }
    void        flush()              { buf.erase( buf.begin(), buf.end()-1); }
    bool        is_space() const { 
	CC_Assert( buf.size() > 0);
	for ( const_iterator i = buf.begin(); i != buf.end() - 1; ++i)
	    if ( *i > ' ')
		return false;
	return true;
    }
};

typedef std::list< Buffer*>         Buffer_list;
typedef Buffer_list::iterator       Buffer_iterator;
typedef Buffer_list::const_iterator Buffer_const_iterator;

void delete_list( Buffer_list* l);

#endif // MODULE_BUFFER //


