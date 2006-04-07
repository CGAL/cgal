// mstring.h
// --------------------------------------------------

// Old mstring.h in mstring_old.h
// We assume a standard conform string class and include it here.
// If it is not available, one may try the mstring_old.h again.

#ifndef MY_STRING_H
#define MY_STRING_H 1

#include <string>

// comment this for SGI 7.2
using namespace std;

// provide hash function for string class
#ifdef __GNUC__

#if ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

#include <ext/hash_map>
#include <ext/hash_set>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_set;

namespace __gnu_cxx {

#else // ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

#include <hash_map>
#include <hash_set>
namespace std {

#endif // ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

template <>
struct hash<string>
{
  size_t operator()(const string& str) const
  {
    unsigned long h = 0;
    const char* s = str.data();
    for (size_t len = str.length(); len > 0; --len, ++s)
      h = 5*h + (unsigned long)(*s);
    return size_t(h);
  }
};

}

#endif // __GNUC__

struct Case_insensitive_string_greater_than {
  int
  compare( const string& s1, const string& s2 ) {
    string::const_iterator it1 = s1.begin();
    string::const_iterator it2 = s2.begin();

    while ( (it1!=s1.end()) && (it2!=s2.end()) )
    {
      if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
      // return -1 to indicate smaller than, 1 otherwise
        return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
      //proceed to the next character in each string
      ++it1;
      ++it2;
    }
    size_t size1=s1.size(), size2=s2.size();// cache lengths
    //return -1,0 or 1 according to strings' lengths
    if (size1==size2)
      return 0;
    return (size1<size2) ? -1 : 1;
  }

  bool
  operator()( const string& s1, const string& s2 ) {
    return compare(s1,s2) > 0;
  }
};

#endif // MY_STRING_H //

// EOF //
