// 
// Smaller changes applied by Lutz Kettner: 
// mainly to conform to the GNU G++ 2.8 compiler and STL.
// Changes are bracketet with comments:
// Begin change L.K.
// End change L.K.
//

// GNU Patch for friend functions (L.K.):
#ifndef GNU_TPRM
#ifdef __GNUC__
#define GNU_TPRM <>
#else
#define GNU_TPRM
#endif
#endif

/**
 **  Copyright (c) 1994-1995 Modena Software Inc.,
 **
 **  Permission to use, copy, modify, distribute and sell this software
 **  and its documentation for any purpose is hereby granted without fee,
 **  provided that the above copyright notice appear in all copies and
 **  that both that copyright notice and this permission notice appear
 **  in supporting documentation.  Modena Software, Inc. makes no 
 **  representations about the suitability of this software for any
 **  purpose.  It is provided "as is" without express or implied warranty.
 **
 **  Copyright (c) 1996 Silicon Graphics
 **
 **  Permission to use, copy, modify, distribute and sell this software
 **  and its documentation for any purpose is hereby granted without fee,
 **  provided that the above copyright notice appear in all copies and
 **  that both that copyright notice and this permission notice appear
 **  in supporting documentation.  Silicon Graphics, Inc. makes no
 **  representations about the suitability of this software for any
 **  purpose.  It is provided "as is" without express or implied warranty.
 **
 **
 **  Changes from the Modena version:
 **
 **  Reference count handling was changed 
 **  to provide for more complete thread-safety.  We now
 **  assume that it is impossible for a reference count to ever overflow
 **  a size_t.  Hence some checks were removed.
 **  
 **  Memory allocation was changed to alloc-style for better performance.
 **
 **  Const_iterator and some other container primitives were added to
 **  improve interaction with STL and to come closer to standard conformance.
 **
 **  Added struct hash<string>.  This requires including hashtable.h,
 **  since hash<string> is a template specialization.
 **
 **  Unchanged from Modena version:
 ** 
 **  It still does not scale well for very long strings.  Building up a million
 **  character string by repeatedly concatenating a character to the end
 **  is likely to take hours on an Indy.
 **/

#ifndef  __cplusplus
#error Must use C++ for BSTRING.H
#endif

#ifndef __MBSTRING_H
#define __MBSTRING_H

#include <ctype.h>
#include <string.h>
#include <iostream.h>
#include <hashtable.h>

#ifdef _SGITHREADS
#include <mutex.h>
#endif

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma set woff 1209
#endif


// bndchk.h
#ifdef BOUNDS_CHECK
void check_bounds
	( 	int index, 
		int container_size, 
		int lineno,
		char *filename	);
#endif

// mexcept.h

#define _THROW_NONE
#define _THROW_DOMAIN
#define _THROW_INVALIDARG
#define _THROW_LENGTH
#define _THROW_OUTRANGE
#define _THROW_RANGE
#define _THROW_OVERFLOW
#define _THROW_ALLOC
#define _THROW_CAST
#define _THROW_TYPEID
#define _THROW_ALLOC_LENGTH
#define _THROW_ALLOC_OUTRANGE
#define _THROW_LENGTH_OUTRANGE
#define _THROW_ALLOC_LENGTH_OUTRANGE


const size_t NPOS  = (size_t)(-1);

enum capacity { default_size, reserve };

template<class charT>
struct  string_char_baggage {

    typedef charT char_type;

    //
    // for users to acquire the basic character type
    //
    // constraints functions
    //
    static void
    assign (char_type& c1, const char_type& c2) _THROW_NONE
    {
        c1 = c2;
    }
    static bool
    eq (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 == c2);
    }
    static bool
    ne (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return !(c1 == c2);
    }
    static bool
    lt (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 < c2);
    }
    static char_type
    eos () _THROW_NONE
    {
        return char_type();     // the null character
    }
    static istream&
    char_in (istream& is, char_type& c) _THROW_NONE
    {
        return is >> c;        // extractor for a char_type object
    }
    static ostream&
    char_out (ostream& os, char_type c) _THROW_NONE
    {
        return os << c;        // inserter for a char_type object
    }
    static bool
    is_del (char_type c) _THROW_NONE
    {
        // characteristic function for delimiters of char_type
        return isspace(c);
    }

    //
    // speed-up functions
    //
    static int
    compare (const char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        for (size_t i = 0; i < n; ++i, ++s1, ++s2)
            if (ne(*s1, *s2))
            {
                return lt(*s1, *s2) ? -1 : 1;
            }
        return 0;
    }
    static size_t
    length (const char_type* s) _THROW_NONE
    {
        size_t l = 0;
        while (ne(*s++, eos()))
            ++l;
        return l;
    }
    static char_type*
    copy (char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        char_type* s = s1;
        for (size_t i = 0; i < n; ++i)
            assign(*++s1, *++s2);
        return s;
    }
};

struct string_char_baggage<char> {

    typedef char char_type;

    //
    // constraint member functions
    //
    static void
    assign (char_type& c1, const char_type& c2) _THROW_NONE
    {
        c1 = c2;
    }
    static bool
    eq (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 == c2);
    }
    static bool
    ne (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 != c2);
    }
    static bool
    lt (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 < c2);
    }
    static char_type
    eos () _THROW_NONE
    {
        return 0;     // the null character
    }
    static istream&
    char_in (istream& is, char_type& c) _THROW_NONE
    {
       // extractor for a char_type object
       // return is >> c;        // this does not work
       is.get(c);
       return is;
    }
    static ostream&
    char_out (ostream& os, char_type c) _THROW_NONE
    {
        return os << c;        // inserter for a char_type object
    }
    static bool
    is_del (char_type c) _THROW_NONE
    {
        // characteristic function for delimiters of char_type
        return isspace(c);
    }

    //
    // speed-up functions
    //
    static int
    compare (const char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        return memcmp(s1, s2, n);
    }
    static size_t
    length (const char_type* s) _THROW_NONE
    {
        return strlen(s);
    }
    static char_type*
    copy (char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        // type cast required as memcpy returns void*
        return (char_type*)memcpy(s1, s2, n);
    }
};

/*
struct string_char_baggage<wchar_t> {

    typedef wchar_t char_type;

    static void
    assign (char_type& c1, const char_type& c2) _THROW_NONE
    {
        c1 = c2;
    }
    static bool
    eq (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 == c2);
    }
    static bool
    ne (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 != c2);
    }
    static bool
    lt (const char_type& c1, const char_type& c2) _THROW_NONE
    {
        return (c1 < c2);
    }
    static char_type
    eos () _THROW_NONE
    {
        return 0;     // the null character
    }
    static istream&
    char_in (istream& is, char_type& c) _THROW_NONE
    {
        return is >> c;        // extractor for a char_type object
    }
    static ostream&
    char_out (ostream& os, char_type c) _THROW_NONE
    {
        return os << c;        // inserter for a char_type object
    }
    static bool
    is_del (char_type c) _THROW_NONE
    {
        // characteristic function for delimiters of char_type
        // using templatized locale::isspace function
        return isspace(c);
    }

    //
    // speed-up functions
    //
    static int
    compare (const char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        return wmemcmp(s1, s2, n);
    }
    static size_t
    length (const char_type* s) _THROW_NONE
    {
        return wcslen(s);
        // May use Koshida's overloaded MSE function strlen(s)
    }
    static char_type*
    copy (char_type* s1, const char_type* s2, size_t n) _THROW_NONE
    {
        return (char_type*)wmemcpy(s1, s2, n);
    }
};
*/

template <class Baggage>
struct __baggage_eq
  : public binary_function<typename Baggage::char_type,
                           typename Baggage::char_type,
                           bool>
{
  bool operator()(const typename Baggage::char_type& c1,
                  const typename Baggage::char_type& c2) const
  {
    return Baggage::eq(c1, c2);
  }
};


template <class charT>
class basic_string;

template <class charT>
class basic_string_ref {

//
// friend class declaration
//
friend class basic_string<charT>;

//
// typedef declarations
//
typedef  string_char_baggage<charT>  baggage_type;

typedef simple_alloc<charT, __ALLOC> charT_alloc;

//
// private data members
//
    charT*   ptr;
    size_t   len;
    size_t   res;
    size_t   count;

// Reference count operations

# ifdef _SGITHREADS
#   if __mips < 3 || !(defined (_ABIN32) || defined(_ABI64))
#   	define __add_and_fetch(l,v) add_then_test((unsigned long *)l,v)
#   endif
    void incr_count ()
    {
	__add_and_fetch(&count, 1);
    }
    size_t decr_count ()
    {
	return __add_and_fetch(&count, (size_t)(-1));
    }
#  else
    void incr_count ()
    {
	++count;
    }
    size_t decr_count ()
    {
	--count;
	return count;
    }
# endif

//
// private constructors and destructors
//
    basic_string_ref () _THROW_NONE ;

    basic_string_ref (size_t size, ::capacity cap) _THROW_ALLOC_LENGTH ;

    basic_string_ref (const basic_string<charT>& str, size_t pos , size_t rlen)
                      _THROW_ALLOC ;

    basic_string_ref (const charT* s, size_t rlen, size_t rres) _THROW_ALLOC ;

    basic_string_ref (const charT* s, size_t n) _THROW_ALLOC_LENGTH ;

    basic_string_ref (const charT* s) _THROW_ALLOC ;

    basic_string_ref (charT c, size_t rep) _THROW_ALLOC_LENGTH ;

    basic_string_ref (const vector<charT>& vec) _THROW_ALLOC_LENGTH ;

    ~basic_string_ref () _THROW_NONE ;

    inline void
    delete_ptr () _THROW_NONE ;

    inline static
    charT
    eos () _THROW_NONE ;

    inline static
    void
    throwlength () _THROW_LENGTH;

    inline static
    void
    throwrange () _THROW_OUTRANGE;

//
//  Replace operator new and delete
//
    void * operator new(size_t n) { return __ALLOC::allocate(n); }

    void operator delete(void *p)
	{ __ALLOC::deallocate(p, sizeof (basic_string_ref<charT>)); }

};

template<class charT>
class basic_string {

private:

//
// typedef declaration
//
    typedef basic_string_ref<charT>   reference_class;
    typedef basic_string_ref<charT>*  reference_pointer;

//
// data member
//
    charT*             c_str_ptr;
    reference_pointer  reference_p;

//
// private member functions
//
    inline charT*
    point () _THROW_NONE ;

    inline size_t&
    len () _THROW_NONE ;

    inline size_t
    ref_count () const _THROW_NONE ;
    
    inline static
    charT
    eos () _THROW_NONE ;

    void
    assign_str (const charT* s, size_t slen) _THROW_ALLOC_LENGTH ;

    void
    append_str (const charT* s, size_t slen) _THROW_ALLOC_LENGTH ;

    void
    insert_str (size_t pos, const charT* s, size_t slen) 
                _THROW_ALLOC_LENGTH_OUTRANGE ;

    void
    replace_str (size_t xlen, size_t pos, const charT* s, size_t slen) 
                 _THROW_ALLOC_LENGTH_OUTRANGE ;

    int
    compare_str (size_t pos, const charT* str, size_t slen, size_t strlen) 
                 const _THROW_OUTRANGE ;

    size_t
    find_str (const charT* s, size_t pos, size_t len) const _THROW_NONE ;

    size_t
    rfind_str (const charT* s, size_t pos, size_t len) const _THROW_NONE ;

    size_t
    find_first_of_str (const charT* s, size_t pos, size_t len) const 
                      _THROW_NONE ;

    size_t
    find_last_of_str (const charT* s, size_t pos, size_t len) const 
                     _THROW_NONE ;

    size_t
    find_first_not_of_str (const charT* s, size_t pos, size_t len) const 
                          _THROW_NONE ;

    size_t
    find_last_not_of_str (const charT* s, size_t pos, size_t len) const 
                         _THROW_NONE ;


protected:

    basic_string (const charT* s, size_t rlen, size_t xlen) _THROW_ALLOC_LENGTH;

    inline void
    delete_ref () _THROW_NONE ;

public:

    typedef  charT                       char_type;
    typedef  string_char_baggage<charT>  baggage_type;


    // Container types
    typedef charT value_type;
    typedef value_type* pointer;
    typedef const value_type* const_iterator;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
// Begin change L.K.
#ifdef __GNUC__
    typedef reverse_iterator<const_iterator>       const_reverse_iterator;
#else
    typedef reverse_iterator<const_iterator, value_type, const_reference,
				 difference_type>  const_reverse_iterator;
#endif
// End change L.K.
    
    // static const size_type npos = NPOS;
    static const size_t npos;

    // iterator support
    const_iterator begin() const { return data(); }
    const_iterator end() const { return data() + reference_p -> len; }
    const_reverse_iterator rbegin() const {
	return const_reverse_iterator(end());
    }
    const_reverse_iterator rend() const {
	return const_reverse_iterator(begin());
    }

    // container support
    size_type size() const { return reference_p -> len; }
    size_type max_size() const { return LONG_MAX; }
    size_type capacity() const
		{ return size_type(reference_p -> res); }
    bool empty() const { return begin() == end(); }

    // constructors
    basic_string () _THROW_ALLOC ;

    basic_string (size_t size, ::capacity cap) _THROW_ALLOC_LENGTH ;

    basic_string (const basic_string<charT>& str, size_t pos = 0, size_t n = NPOS)
                  _THROW_ALLOC_OUTRANGE ;

    basic_string (const charT* s, size_t n) _THROW_ALLOC_LENGTH ;

    basic_string (const charT* s) _THROW_ALLOC ;

    basic_string (charT c, size_t rep = 1) _THROW_ALLOC_LENGTH ;

    basic_string (const vector<charT>& vec) _THROW_ALLOC_LENGTH ;

    ~basic_string () _THROW_NONE ;

    basic_string<charT>&
    operator= (const basic_string<charT>& str) _THROW_ALLOC ;

    basic_string<charT>&
    operator= (const charT* s) _THROW_ALLOC ;

    basic_string<charT>&
    operator= (charT c) _THROW_ALLOC ;

    basic_string<charT>&
    operator+= (const basic_string<charT>& rhs) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    operator+= (const charT* s) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    operator+= (charT c) _THROW_ALLOC_LENGTH ;

    operator vector<charT> () const _THROW_ALLOC { 
        return vector<charT> (data(), data()+length());
	}

    basic_string<charT>&
    append (const basic_string<charT>& str, size_t pos = 0, size_t n = NPOS)
            _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    append (const charT* s, size_t n) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    append (const charT* s) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    append (charT c, size_t rep = 1) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    assign (const basic_string<charT>& str, size_t pos = 0, size_t n = NPOS)
            _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    assign (const charT* s, size_t n) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    assign (const charT* s) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    assign (charT c, size_t rep = 1) _THROW_ALLOC_LENGTH ;

    basic_string<charT>&
    insert (size_t pos1, const basic_string<charT>& str, size_t pos2 = 0,
            size_t n = NPOS) _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    insert (size_t pos, const charT* s, size_t n) _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    insert (size_t pos, const charT* s) _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    insert (size_t pos, charT c, size_t rep = 1) _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    remove (size_t pos = 0, size_t n = NPOS) _THROW_ALLOC_OUTRANGE ;

    basic_string<charT>&
    replace (size_t pos1, size_t n1, const basic_string<charT>& str, size_t pos2 = 0,
             size_t n2 = NPOS) _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    replace (size_t pos, size_t n1, const charT* s, size_t n2)
             _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    replace (size_t pos, size_t n1, const charT* s)
             _THROW_ALLOC_LENGTH_OUTRANGE ;

    basic_string<charT>&
    replace (size_t pos, size_t n, charT c, size_t rep = 1)
             _THROW_ALLOC_LENGTH_OUTRANGE ;

    inline charT
    at (size_t pos) const _THROW_OUTRANGE ;

    void
    put_at (size_t pos, charT c) _THROW_ALLOC_OUTRANGE ;

    inline charT
    operator[] (size_t pos) const _THROW_NONE ;

    charT&
    operator[] (size_t pos) _THROW_ALLOC_OUTRANGE ;

    const charT*
    c_str () const _THROW_ALLOC ;

    inline const charT*
    data () const _THROW_NONE ;

    inline size_t
    length () const _THROW_NONE ;

    void
    resize (size_t n, charT c) _THROW_ALLOC_LENGTH ;

    void
    resize (size_t n) _THROW_ALLOC_LENGTH ;

    inline size_t
    reserve () const _THROW_NONE ;

    void
    reserve (size_t res_arg) _THROW_ALLOC_LENGTH ;

    size_t
    copy (charT* s, size_t n, size_t pos = 0) const _THROW_OUTRANGE ;

    size_t
    find (const basic_string<charT>& str, size_t pos = 0) const _THROW_NONE ;

    size_t
    find (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    find (const charT* s, size_t pos = 0) const _THROW_NONE ;

    size_t
    find (charT c, size_t pos = 0) const _THROW_NONE ;

    size_t
    rfind (const basic_string<charT>& str, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    rfind (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    rfind (const charT* s, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    rfind (charT c, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    find_first_of (const basic_string<charT>& str, size_t pos = 0) const _THROW_NONE ;

    size_t
    find_first_of (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    find_first_of (const charT* s, size_t pos = 0) const _THROW_NONE ;

    size_t
    find_first_of (charT c, size_t pos = 0) const _THROW_NONE ;


    size_t
    find_last_of (const basic_string<charT>& str, size_t pos = NPOS) const
                  _THROW_NONE ;

    size_t
    find_last_of (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    find_last_of (const charT* s, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    find_last_of (charT c, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    find_first_not_of (const basic_string<charT>& str, size_t pos = 0) const
                       _THROW_NONE ;

    size_t
    find_first_not_of (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    find_first_not_of (const charT* s, size_t pos = 0) const _THROW_NONE ;

    size_t
    find_first_not_of (charT c, size_t pos = 0) const _THROW_NONE ;

    size_t
    find_last_not_of (const basic_string<charT>& str, size_t pos = NPOS) const
                      _THROW_NONE ;

    size_t
    find_last_not_of (const charT* s, size_t pos, size_t n) const _THROW_NONE ;

    size_t
    find_last_not_of (const charT* s, size_t pos = NPOS) const _THROW_NONE ;

    size_t
    find_last_not_of (charT c, size_t pos = NPOS) const _THROW_NONE ;

    basic_string<charT>
    substr (size_t pos = 0,  size_t n = NPOS) const _THROW_ALLOC_OUTRANGE ;

    int
    compare (const basic_string<charT>& str, size_t pos = 0, size_t n = NPOS) const
             _THROW_OUTRANGE ;

    int
    compare (const charT* s, size_t pos, size_t n) const
             _THROW_LENGTH_OUTRANGE ;

    int
    compare (const charT* s, size_t pos = 0) const _THROW_OUTRANGE ;

    int
    compare (charT c, size_t pos = 0, size_t rep = 1) const 
             _THROW_LENGTH_OUTRANGE ;

    friend
    ostream&
    operator<< GNU_TPRM(ostream& o, const basic_string<charT>& s) _THROW_NONE ;

    friend
    istream&
    operator>> GNU_TPRM(istream& i, basic_string<charT>& s)
    _THROW_ALLOC_LENGTH ;

    friend
    basic_string<charT>
    operator+ GNU_TPRM(const basic_string<charT>& lhs, 
		       const basic_string<charT>& rhs)
    _THROW_ALLOC_LENGTH ;

    friend
    basic_string<charT>
    operator+ GNU_TPRM(const charT* lhs, const basic_string<charT>& rhs)
    _THROW_ALLOC_LENGTH ;

    friend
    basic_string<charT>
    operator+ GNU_TPRM(charT lhs, const basic_string<charT>& rhs) 
    _THROW_ALLOC_LENGTH ;

    friend
    basic_string<charT>
    operator+ GNU_TPRM(const basic_string<charT>& lhs, const charT* rhs)
    _THROW_ALLOC_LENGTH ;

    friend
    basic_string<charT>
    operator+ GNU_TPRM(const basic_string<charT>& lhs, charT rhs) 
    _THROW_ALLOC_LENGTH ;

};

template<class charT>
const size_t basic_string<charT>::npos = NPOS;

template <class charT>
inline void
basic_string_ref<charT>::delete_ptr () _THROW_NONE
{
    if (res)
    {
	charT_alloc::deallocate(ptr,res);
        res = 0;
        ptr = 0;
    }
}


template <class charT>
inline void
basic_string_ref<charT>::throwlength ()  _THROW_LENGTH
{
//#ifdef  __MEXCEPT
//   throw   length_error("Length exception occurred");
//#else
   cout << "Length exception occurred" << endl;
   abort();
//#endif
}

template <class charT>
inline void
basic_string_ref<charT>::throwrange () _THROW_OUTRANGE
{
//#ifdef  __MEXCEPT
//   throw   out_of_range("Out of range exception occurred");
//#else
   cout << "Out of range exception occurred" << endl;
   abort();
//#endif
}


template <class charT>
inline void
basic_string<charT>::delete_ref () _THROW_NONE
{
    if (!(reference_p->decr_count()))
       delete reference_p;
}

template <class charT>
inline size_t
basic_string<charT>::ref_count () const _THROW_NONE
{
    return reference_p->count;
}

template <class charT>
inline const charT*
basic_string<charT>::data () const _THROW_NONE
{
    if (length())
        return reference_p->ptr;
    else
        return 0;
}

template <class charT>
inline charT*
basic_string<charT>::point () _THROW_NONE
{
    return reference_p->ptr;
}

template <class charT>
inline size_t&
basic_string<charT>::len () _THROW_NONE
{
    return reference_p->len;
}

template <class charT>
inline size_t
basic_string<charT>::length () const _THROW_NONE
{
    return reference_p->len;
}

template <class charT>
inline size_t
basic_string<charT>::reserve () const _THROW_NONE
{
    return reference_p->res;
}

template <class charT>
inline charT
basic_string<charT>::at (size_t pos) const _THROW_OUTRANGE
{
    if (pos >= length())
    {
        reference_class::throwrange();
    }
    return *(data()+pos);
}

template <class charT>
inline charT
basic_string<charT>::operator[] (size_t pos) const _THROW_NONE
{
    if (pos < length())
        return *(data()+pos);
    else
        return 0;
}

template <class charT>
inline bool
operator== (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
            _THROW_NONE
{
    return !(lhs.compare(rhs));
}

template <class charT>
inline bool
operator== (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(rhs.compare(lhs));
}

template <class charT>
inline bool
operator== (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(rhs.compare(lhs));
}

template <class charT>
inline bool
operator== (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    return !(lhs.compare(rhs));
}

template <class charT>
inline bool
operator== (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    return !(lhs.compare(rhs));
}

#ifdef __MNONDEF
template <class charT>
inline bool
operator!= (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
            _THROW_NONE
{
    return lhs.compare(rhs);
}
#endif

template <class charT>
inline bool
operator!= (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return rhs.compare(lhs);
}

template <class charT>
inline bool
operator!= (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return rhs.compare(lhs);
}

template <class charT>
inline bool
operator!= (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    return lhs.compare(rhs);
}

template <class charT>
inline bool
operator!= (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    return lhs.compare(rhs);
}

template <class charT>
inline bool
operator< (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
           _THROW_NONE
{
    if (lhs.compare(rhs) < 0)
        return true;
    else
        return false;
}

template <class charT>
inline bool
operator< (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    if (rhs.compare(lhs) > 0)
        return true;
    else
        return false;
}

template <class charT>
inline bool
operator< (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    if (rhs.compare(lhs) > 0)
        return true;
    else
        return false;
}

template <class charT>
inline bool
operator< (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    if (lhs.compare(rhs) < 0)
        return true;
    else
        return false;
}

template <class charT>
inline bool
operator< (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    if (lhs.compare(rhs) < 0)
        return true;
    else
        return false;
}

#ifdef  __MNONDEF
template <class charT>
inline bool
operator> (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
           _THROW_NONE
{
    return (rhs < lhs);
}
#endif

template <class charT>
inline bool
operator> (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return (rhs < lhs);
}

template <class charT>
inline bool
operator> (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return (rhs < lhs);
}

template <class charT>
inline bool
operator> (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    return (rhs < lhs);
}

template <class charT>
inline bool
operator> (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    return (rhs < lhs);
}

#ifdef  __MNONDEF
template <class charT>
inline bool
operator>= (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
            _THROW_NONE
{
    return !(lhs < rhs);
}
#endif

template <class charT>
inline bool
operator>= (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(lhs < rhs);
}

template <class charT>
inline bool
operator>= (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(lhs < rhs);
}

template <class charT>
inline bool
operator>= (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    return !(lhs < rhs);
}

template <class charT>
inline bool
operator>= (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    return !(lhs < rhs);
}

#ifdef  __MNONDEF
template <class charT>
inline bool
operator<= (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
            _THROW_NONE
{
    return !(rhs < lhs);
}
#endif

template <class charT>
inline bool
operator<= (const charT* lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(rhs < lhs);
}

template <class charT>
inline bool
operator<= (charT lhs, const basic_string<charT>& rhs) _THROW_NONE
{
    return !(rhs < lhs);
}

template <class charT>
inline bool
operator<= (const basic_string<charT>& lhs, const charT* rhs) _THROW_NONE
{
    return !(rhs < lhs);
}

template <class charT>
inline bool
operator<= (const basic_string<charT>& lhs, charT rhs) _THROW_NONE
{
    return !(rhs < lhs);
}

// definitions : can be in a .c file
//

template <class charT>
charT
basic_string_ref<charT>::eos () _THROW_NONE
{
    return baggage_type::eos();
}

template <class charT>
basic_string_ref<charT>::basic_string_ref () _THROW_NONE
{
     res = len = 0;
     ptr = 0;
     count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (size_t size, ::capacity cap)
                                           _THROW_ALLOC_LENGTH
{
    if (cap == ::reserve)
    {
       len = 0;
       res = size;
       ptr = charT_alloc::allocate(res);
    }
    else if ((cap == ::default_size) && (size != NPOS))
    {
       res = len = size;
       if (res)
       {
          ptr = charT_alloc::allocate(res);
          for (size_t position = 0; position < len; ++position)
              baggage_type::assign (*(ptr+position), eos());
       }
       else
          ptr = 0;
    }
    else
    {
       throwlength();
    }
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (const basic_string<charT>& str,
                                           size_t pos, size_t rlen) _THROW_ALLOC
{
    res = len = rlen;
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       baggage_type::copy (ptr, str.data()+pos, len);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (const charT* s, size_t rlen,
                                           size_t rres) _THROW_ALLOC
{
    res = rres;
    len = rlen;
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       if (len)
          baggage_type::copy (ptr, s, len);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (const charT* s, size_t n)
                                           _THROW_ALLOC_LENGTH
{
    if (n == NPOS)
    {
       throwlength();
    }
    res = len = n;
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       baggage_type::copy (ptr, s, len);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (const charT* s) _THROW_ALLOC
{
    res = len = baggage_type::length(s);
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       baggage_type::copy (ptr, s, len);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (charT c, size_t rep)
                                           _THROW_ALLOC_LENGTH
{
    if (rep == NPOS)
    {
       throwlength();
    }
    res = len = rep;
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       for (size_t  position = 0; position < len; ++position)
            baggage_type::assign (*(ptr+position), c);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::basic_string_ref (const vector<charT>& vec)
                                           _THROW_ALLOC_LENGTH
{
    size_t  n = vec.size();
    if (n == NPOS)
    {
        throwlength();
    }
    res = len = n;
    if (res)
    {
       ptr = charT_alloc::allocate(res);
       baggage_type::copy (ptr, vec.begin(), len);
    }
    else
       ptr = 0;
    count = 1;
}

template <class charT>
basic_string_ref<charT>::~basic_string_ref () _THROW_NONE
{
    delete_ptr();
}

template <class charT>
charT
basic_string<charT>::eos () _THROW_NONE
{
    return baggage_type::eos();
}

template <class charT>
void
basic_string<charT>::assign_str (const charT* s, size_t slen)
                                 _THROW_ALLOC_LENGTH
{
    if (slen == NPOS)
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (slen && (reserve() < slen)))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (s, slen);
        delete_ref();
        reference_p = tmp;
    }
    else if (slen)
    {
        baggage_type::copy (point(), s, slen);
    }
    reference_p->len = slen;
}

template <class charT>
void
basic_string<charT>::append_str (const charT* s, size_t slen)
                                 _THROW_ALLOC_LENGTH
{
    if (length() >= (NPOS-slen))
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (slen > (reserve()-length())))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), length(), length()+slen);
        delete_ref();
        reference_p = tmp;
    }
    if (slen)
        baggage_type::copy (point()+length(), s, slen);
    reference_p->len += slen;
}

template <class charT>
void
basic_string<charT>::insert_str (size_t pos, const charT* s, size_t slen)
                                 _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if (length() >= (NPOS-slen))
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (slen > (reserve()-length())))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), pos, length()+slen);
        baggage_type::copy (tmp->ptr+pos+slen, data()+pos, length()-pos);
        tmp->len = length();
        delete_ref();
        reference_p = tmp;
    }
    else
    {
        for (size_t count = length()-pos; count > 0; --count)
             baggage_type::assign (*(point()+pos+slen+count-1),
                                   *(data()+pos+count-1));
    }
    if (slen)
        baggage_type::copy (point()+pos, s, slen);
    reference_p->len += slen;
}

template <class charT>
void
basic_string<charT>::replace_str (size_t xlen, size_t pos, const charT* s,
                     size_t slen) _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if ((length()-xlen) >= (NPOS-slen))
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (reserve() < (length()+slen-xlen)))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), pos, length()+slen-xlen);
        baggage_type::copy (tmp->ptr+pos+slen, data()+pos+xlen,
                            length()-pos-xlen);
        tmp->len = length();
        delete_ref();
        reference_p = tmp;
    }
    else 
    {
        if (slen < xlen)
            baggage_type::copy (point()+pos+slen, data()+pos+xlen,
                                length()-pos-xlen);
        else
        {
            for (size_t count = length()-pos-xlen; count > 0; --count)
                baggage_type::assign (*(point()+pos+slen+count-1),
                                      *(data()+pos+xlen+count-1));
        }
    }
    if (slen)
        baggage_type::copy (point()+pos, s, slen);
    reference_p->len += (slen-xlen);
}

template <class charT>
int
basic_string<charT>::compare_str (size_t pos, const charT* str, size_t slen,
                     size_t strlen) const _THROW_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    size_t rlen = (slen > strlen ) ?  strlen  : slen;
    int result;
    if (!length())
        return str ? (eos()- *str) : eos();
    result =  baggage_type::compare (data()+pos, str, rlen);
    return result ? result : (length()-pos-strlen);
}

template <class charT>
size_t
basic_string<charT>::find_str (const charT* s, size_t pos, size_t len) 
                               const _THROW_NONE
{
    size_t  count = pos;
    size_t  shift;
    size_t  place;
    if ((length() == 0) || (len == 0))
        return NPOS;
    while (len <= (length()-count))
    {
        for (place = 0; place < len; ++place)
        {
            if (baggage_type::ne(*(s+len-1-place),
				 *(data()+count+(len-1-place))))
                break;
        }
        if (place == len)
            return count;
        shift = find(*(s+len-1-place), count+(len-place));
        if (shift == NPOS)
            return NPOS;
        count = shift-(len-place-1);
    }
    return NPOS;
}

template <class charT>
size_t
basic_string<charT>::rfind_str (const charT* s, size_t pos, size_t len) 
                                const _THROW_NONE
{
// Begin change L.K.
#ifdef __GNUC__
    typedef reverse_iterator<const charT*>       const_reverse_ptr;
#else
    typedef reverse_iterator<const charT*, charT, const charT&, ptrdiff_t>
            const_reverse_ptr;
#endif
// End change L.K.

    size_t N = length();
    if (N < len || len == 0)
      return NPOS;
    else {
      const_iterator last = begin() + ((pos >= N - len) ? N : pos + len);
      const_reverse_iterator rfirst = const_reverse_iterator(last);
      const_reverse_iterator rlast = const_reverse_iterator(begin());

      const_reverse_iterator rit =
        search(rfirst, rlast,
               const_reverse_ptr(s + len), const_reverse_ptr(s),
               __baggage_eq<baggage_type>());
      if (rit == rlast)
        return NPOS;
      else
        return (rit.base() - len) - begin();      
    }
}


template <class charT>
size_t
basic_string<charT>::find_first_of_str (const charT* s, size_t pos, size_t len) 
                                        const _THROW_NONE
{
    size_t temp;
    size_t count  = pos;
    size_t result = NPOS;
    while (count < length())
    {
        temp = 0;
        while ((temp < len) && baggage_type::ne(*(data()+count), *(s+temp)))
            ++temp;
        if (temp != len)
            break;
        ++count;
    }
    temp = (count >= length()) ? NPOS : count;
    return ((result > temp) ? temp : result);
}

template <class charT>
size_t
basic_string<charT>::find_last_of_str (const charT* s, size_t pos, size_t len) 
                                       const _THROW_NONE
{
    size_t temp = 0;
    size_t count = (pos < length()) ? (pos+1) : length();
    if (length())
    {
       while (count > 0)
       {
           temp = 0;
           --count;
           while ((temp != len) && baggage_type::ne(*(data()+count), *(s+temp)))
               ++temp;
           if (temp != len)
               break;
       }
    }
    return ((temp != len) && length()) ? count : NPOS;
}

template <class charT>
size_t
basic_string<charT>::find_first_not_of_str (const charT* s, size_t pos,
                     size_t len) const _THROW_NONE
{
    size_t count = pos;
    while (count < length())
    {
        size_t temp  = 0;
        while (temp < len)
        {
            if (baggage_type::eq(*(data()+count), *(s+temp)))
                break;
            ++temp;
        }
        if (temp == len)
            break;
        ++count;
    }
    return  ((count >= length()) ? NPOS : count);
}

template <class charT>
size_t
basic_string<charT>::find_last_not_of_str (const charT* s, size_t pos,
                     size_t len) const _THROW_NONE
{
    size_t temp = 0;
    size_t count = (pos < length()) ? (pos+1) : length();

    if (length())
    {
        while (count > 0)
        {
           temp = 0;
           while (temp != len)
           {
               if (baggage_type::eq(*(data()+count-1), *(s+temp)))
                   break;
               ++temp;
           }
           if (temp == len)
               break;
           --count;
       }
    }
    return ((temp == len) && length()) ? count-1 : NPOS;
}

template <class charT>
basic_string<charT>::basic_string () _THROW_ALLOC
{
    reference_p = new basic_string_ref<charT> ();
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (size_t size, ::capacity cap)
                                   _THROW_ALLOC_LENGTH
{
    reference_p = new basic_string_ref<charT> (size, cap);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (const basic_string<charT>& str,
                                   size_t pos, size_t n) _THROW_ALLOC_OUTRANGE
{
    if (pos > str.length())
    {
       reference_class::throwrange();
    }
    size_t rlen =  (n > (str.length() - pos)) ? str.length() - pos : n;
    if (rlen == str.length())
       (reference_p = str.reference_p)->incr_count();
    else
        reference_p = new basic_string_ref<charT> (str, pos, rlen);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (const charT* s, size_t rlen, size_t xlen)
                                   _THROW_ALLOC_LENGTH
{
    if (rlen >= (NPOS - xlen))
    {
        reference_class::throwlength();
    }
    reference_p = new basic_string_ref<charT> (s, rlen, rlen+xlen);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (const charT* s, size_t n) _THROW_ALLOC_LENGTH
{
    reference_p = new basic_string_ref<charT> (s, n);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (const charT* s) _THROW_ALLOC
{
    reference_p = new basic_string_ref<charT> (s);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (charT c, size_t rep) _THROW_ALLOC_LENGTH
{
    reference_p = new basic_string_ref<charT> (c, rep);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::basic_string (const vector<charT>& vec) _THROW_ALLOC_LENGTH
{
    reference_p = new basic_string_ref<charT> (vec);
    c_str_ptr = 0;
}

template <class charT>
basic_string<charT>::~basic_string () _THROW_NONE
{
    delete_ref();
    if (c_str_ptr)
        delete[] c_str_ptr;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator= (const basic_string<charT>& str) _THROW_ALLOC
{
    if (this != &str)
    {
        delete_ref();
        (reference_p = str.reference_p)->incr_count();
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator= (const charT* s) _THROW_ALLOC
{
    assign_str (s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator= (charT c) _THROW_ALLOC
{
    if ((ref_count() == 1) && (reserve() >= 1))
    {
        baggage_type::assign (*(point()), c);
        reference_p->len = 1;
    }
    else
    {
        delete_ref();
        reference_p = new basic_string_ref<charT> (c, 1);
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator+= (const basic_string<charT>& rhs) _THROW_ALLOC_LENGTH
{
    append_str (rhs.data(), rhs.length());
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator+= (const charT* s) _THROW_ALLOC_LENGTH
{
    append_str (s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::operator+= (charT c) _THROW_ALLOC_LENGTH
{
    if (length() >= (NPOS-1))
    {
        reference_class::throwlength();
    }
    if (!((ref_count() == 1) && (reserve() > length())))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), length(), length()+1);
        delete_ref();
        reference_p = tmp;
    }
    baggage_type::assign (*(point()+length()), c);
    reference_p->len++;
    return *this;
}


template <class charT>
basic_string<charT>&
basic_string<charT>::append (const basic_string<charT>& str, size_t pos, size_t n)
                     _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > str.length())
    {
        reference_class::throwrange();
    }
    append_str (str.data() + pos, (n>(str.length()-pos))?(str.length()-pos):n);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::append (const charT* s, size_t n) _THROW_ALLOC_LENGTH
{
    append_str (s, n);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::append (const charT* s) _THROW_ALLOC_LENGTH
{
    append_str (s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::append (charT c, size_t rep) _THROW_ALLOC_LENGTH
{
    if (length() >= (NPOS-rep))
    {
        reference_class::throwlength();
    }
    if (rep)
    {
       if ((ref_count() > 1) || (reserve() < (length() + rep)))
       {
          reference_pointer tmp;
          tmp = new basic_string_ref<charT> (data(), length(), length()+rep);
          delete_ref();
          reference_p = tmp;
       }
       for (size_t count = 0; count < rep; ++count)
            baggage_type::assign (*(point()+length()+count), c);
       reference_p->len += rep;
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::assign (const basic_string<charT>& str, size_t pos, size_t n)
                             _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > str.length())
    {
        reference_class::throwrange();
    }
    size_t rlen = (n > (str.length() - pos)) ? str.length() - pos : n;
    if ((rlen == str.length()) && (str.ref_count() != NPOS))
    {
       delete_ref();
       (reference_p = str.reference_p)->count++;
    }
    else
       assign_str (str.data()+pos, rlen);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::assign (const charT* s, size_t n) _THROW_ALLOC_LENGTH
{
    assign_str (s, n);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::assign (const charT* s) _THROW_ALLOC_LENGTH
{
    assign_str (s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::assign (charT c, size_t rep) _THROW_ALLOC_LENGTH
{
    if (rep == NPOS)
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (rep && (reserve() < rep)))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (c, rep);
        delete_ref();
        reference_p = tmp;
    }
    else
    {
        for (size_t count = 0; count < rep; ++count)
            baggage_type::assign (*(point()+count), c);
        reference_p->len = rep;
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::insert (size_t pos1, const basic_string<charT>& str,
                             size_t pos2, size_t n) _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos2 > str.length())
    {
        reference_class::throwrange();
    }
    size_t rlen = (n > (str.length() - pos2)) ? str.length() - pos2 : n;
    insert_str (pos1, str.data()+pos2, rlen);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::insert (size_t pos, const charT* s, size_t n)
                             _THROW_ALLOC_LENGTH_OUTRANGE
{
    insert_str(pos, s, n);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::insert (size_t pos, const charT* s)
                             _THROW_ALLOC_LENGTH_OUTRANGE
{
    insert_str(pos, s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::insert (size_t pos, charT c, size_t rep)
                             _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if ((rep == NPOS) || (length() >= (NPOS - rep)))
    {
        reference_class::throwlength();
    }
    if (rep)
    {
        size_t count;
        if ((ref_count() > 1) || (reserve() < (length()+rep)))
        {
           reference_pointer tmp;
           tmp = new basic_string_ref<charT> (data(), pos, length()+rep);
           if (length())
               for (count = length()-pos; count > 0; --count)
                    baggage_type::assign (*(tmp->ptr+pos+rep+count-1),
                                          *(data()+pos+count-1));
           tmp->len = length();
           delete_ref();
           reference_p = tmp;
        }
        else
        {
           for (count = length()-pos; count > 0; --count)
                baggage_type::assign (*(point()+pos+rep+count-1),
                                      *(data()+pos+count-1));
        }
        for (count = 0; count < rep; ++count)
            baggage_type::assign (*(point()+pos+count), c);
        reference_p->len += rep;
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::remove (size_t pos, size_t n) _THROW_ALLOC_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    size_t xlen = (n > (length()-pos)) ? (length()-pos) : n;
    if (ref_count() > 1)
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), pos, length());
        baggage_type::copy (tmp->ptr+pos, data()+pos+xlen, length()-xlen-pos);
        tmp->len = length()-xlen;
        delete_ref();
        reference_p = tmp;
    }
    else if (xlen == length())
        reference_p->len = 0;
    else if (xlen)
    {
        baggage_type::copy (point()+pos, data()+pos+xlen, length()-xlen-pos);
        reference_p->len -= xlen;
    }
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::replace (size_t pos1, size_t n1, const basic_string<charT>& str,
                     size_t pos2, size_t n2) _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos2 > str.length())
    {
        reference_class::throwrange();
    }
    size_t xlen = (n1 > (length()-pos1)) ? (length()-pos1) : n1;
    size_t rlen = (n2 > (str.length()-pos2)) ? (str.length()-pos2) : n2;
    replace_str (xlen, pos1, str.data()+pos2, rlen);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::replace (size_t pos, size_t n1, const charT* s, size_t n2)
                              _THROW_ALLOC_LENGTH_OUTRANGE
{
    size_t xlen = (n1 > (length()-pos)) ? (length()-pos) : n1;
    replace_str (xlen, pos, s, n2);
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::replace (size_t pos, size_t n1, const charT* s)
                              _THROW_ALLOC_LENGTH_OUTRANGE
{
    size_t xlen = (n1 > (length()-pos)) ? (length()-pos) : n1;
    replace_str (xlen, pos, s, baggage_type::length(s));
    return *this;
}

template <class charT>
basic_string<charT>&
basic_string<charT>::replace (size_t pos, size_t n, charT c, size_t rep)
                              _THROW_ALLOC_LENGTH_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    size_t xlen = (n > (length()-pos)) ? (length()-pos) : n;
    if ((length()-xlen) >= (NPOS-rep))
    {
        reference_class::throwlength();
    }
    if (!rep)
        return remove (pos, n);
    else
    {
        size_t count;
        if ((ref_count() > 1) || (reserve() < (length()-xlen+rep)))
        {
            reference_pointer tmp;
            tmp = new basic_string_ref<charT> (data(), pos,
                  length()+((xlen > rep) ? (xlen-rep) : 0));
            if (rep < xlen)
                baggage_type::copy (tmp->ptr+pos+rep, data()+pos+xlen,
                                    length()-pos-xlen);
            else
            {
                for (count = length()-xlen-pos; count > 0; --count)
                    baggage_type::assign (*(tmp->ptr+pos+rep+count-1),
                                          *(data()+pos+xlen+count-1));
            }
            tmp->len = length();
            delete_ref();
            reference_p = tmp;
        }
        else
        {
            if (rep < xlen)
                baggage_type::copy (point()+pos+rep, data()+pos+xlen,
                                    length()-pos-xlen);
            else
            {
                for (count = length()-xlen-pos; count > 0; --count)
                    baggage_type::assign (*(point()+pos+rep+count-1),
                                          *(data()+pos+xlen+count-1));
            }
        }
        for (count = 0; count < rep; ++count)
            baggage_type::assign (*(point()+pos+count), c);
        reference_p->len += (rep-xlen);
    }
    return *this;
}

template <class charT>
void
basic_string<charT>::put_at (size_t pos, charT c) _THROW_ALLOC_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if ((ref_count() > 1) || (pos == reserve()))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), length(),
              length()+((pos==length())?1:0));
        delete_ref();
        reference_p = tmp;
    }
    if (pos == length())
        ++reference_p->len;
    baggage_type::assign (*(point()+pos), c);
}

template <class charT>
charT&
basic_string<charT>::operator[] (size_t pos) _THROW_ALLOC_OUTRANGE
{
    if (pos >= length())
    {
        reference_class::throwrange();
    }
    if (ref_count() > 1)
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(), length(), length());
        delete_ref();
        reference_p = tmp;
    }
    return *(point()+pos);
}

template <class charT>
const charT*
basic_string<charT>::c_str () const _THROW_ALLOC
{
    if (c_str_ptr)
        delete[] ((basic_string<charT>*)this)->c_str_ptr;
    ((basic_string<charT>*)this)->c_str_ptr = new charT [length()+1];
    if (length())
        baggage_type::copy (((basic_string<charT>*)this)->c_str_ptr, data(),
                              length());
    baggage_type::assign (*(((basic_string<charT>*)this)->c_str_ptr+length()),
                          eos());
    return c_str_ptr;
}

template <class charT>
void
basic_string<charT>::resize (size_t n, charT c) _THROW_ALLOC_LENGTH
{
    if (n == NPOS)
    {
        reference_class::throwlength();
    }
    if ((ref_count() > 1) || (n > reserve()))
    {
        reference_pointer tmp;
        tmp = new basic_string_ref<charT> (data(),
              ((n > length()) ? length() : n), n);
        delete_ref();
        reference_p = tmp;
    }
    while (reference_p->len < n)
    {
        baggage_type::assign (*(reference_p->ptr+length()), c);
        ++reference_p->len;
    }
    reference_p->len = n;
}

template <class charT>
void
basic_string<charT>::resize (size_t n) _THROW_ALLOC_LENGTH
{
    resize (n, eos());
}

template <class charT>
void
basic_string<charT>::reserve (size_t res_arg) _THROW_ALLOC_LENGTH
{
    if (res_arg == NPOS)
    {
        reference_class::throwlength();
    }
    if (res_arg > reserve())
    {
        reference_pointer tmp;
        tmp = new  basic_string_ref<charT> (data(), length(), res_arg);
        delete_ref();
        reference_p = tmp;
    }
}

template <class charT>
size_t
basic_string<charT>::copy (charT* s, size_t n, size_t pos) const _THROW_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    size_t  rlen = (n > (length()-pos)) ? (length()-pos) : n;
    if (length())
        baggage_type::copy (s, data()+pos, rlen);
    return rlen;
}

template <class charT>
size_t
basic_string<charT>::find (const basic_string<charT>& str, size_t pos) const
                           _THROW_NONE
{
    return find_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::find (const charT* s, size_t pos, size_t n) const
                           _THROW_NONE
{
    return find_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::find (const charT* s, size_t pos) const _THROW_NONE
{
    return find_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::find (charT c, size_t pos) const _THROW_NONE
{
    while ((pos < length()) && (baggage_type::ne(*(data()+pos), c)))
        ++pos;
    return ((pos < length()) ? pos : NPOS);
}

template <class charT>
size_t
basic_string<charT>::rfind (const basic_string<charT>& str, size_t pos) const
                            _THROW_NONE
{
    return rfind_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::rfind (const charT* s, size_t pos, size_t n) const
                            _THROW_NONE
{
    return rfind_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::rfind (const charT* s, size_t pos) const _THROW_NONE
{
    return rfind_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::rfind (charT c, size_t pos) const _THROW_NONE
{
    size_t count = ((pos < length()) ? pos+1 : length());
    if (length() == 0)
        return NPOS;
    while ((baggage_type::ne(*(data()+count-1), c)) && (count > 1))
        --count;
    if ((count == 1) && (baggage_type::ne(*(data()), c)))
        return NPOS;
    else
        return count-1;
}

template <class charT>
size_t
basic_string<charT>::find_first_of (const basic_string<charT>& str, size_t pos) const
                                    _THROW_NONE
{
    return find_first_of_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::find_first_of (const charT* s, size_t pos, size_t n) const
                                    _THROW_NONE
{
    return find_first_of_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::find_first_of (const charT* s, size_t pos) const
                                    _THROW_NONE
{
    return find_first_of_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::find_first_of (charT c, size_t pos) const _THROW_NONE
{
    return find (c, pos);
}

template <class charT>
size_t
basic_string<charT>::find_last_of (const basic_string<charT>& str, size_t pos) const
                                   _THROW_NONE
{
    return find_last_of_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::find_last_of (const charT* s, size_t pos, size_t n) const
                                   _THROW_NONE
{
    return find_last_of_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::find_last_of (const charT* s, size_t pos) const _THROW_NONE
{
    return find_last_of_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::find_last_of (charT c, size_t pos) const _THROW_NONE
{
    return rfind (c, pos);
}

template <class charT>
size_t
basic_string<charT>::find_first_not_of (const basic_string<charT>& str, size_t pos)
                                        const _THROW_NONE
{
    return find_first_not_of_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::find_first_not_of (const charT* s, size_t pos, size_t n)
                                        const _THROW_NONE
{
    return find_first_not_of_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::find_first_not_of (const charT* s, size_t pos) const
                                        _THROW_NONE
{
    return find_first_not_of_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::find_first_not_of (charT c, size_t pos) const _THROW_NONE
{
    while ((pos < length()) && (baggage_type::eq(*(data()+pos), c)))
        ++pos;
    return ((pos < length()) ? pos : NPOS);
}

template <class charT>
size_t
basic_string<charT>::find_last_not_of (const basic_string<charT>& str, size_t pos)
                                       const _THROW_NONE
{
    return find_last_not_of_str (str.data(), pos, str.length());
}

template <class charT>
size_t
basic_string<charT>::find_last_not_of (const charT* s, size_t pos, size_t n)
                                       const _THROW_NONE
{
    return find_last_not_of_str (s, pos, n);
}

template <class charT>
size_t
basic_string<charT>::find_last_not_of (const charT* s, size_t pos) const
                                       _THROW_NONE
{
    return find_last_not_of_str (s, pos, baggage_type::length(s));
}

template <class charT>
size_t
basic_string<charT>::find_last_not_of (charT c, size_t pos) const _THROW_NONE
{
    size_t count = ((pos < length()) ? pos+1 : length());
    if (length() == 0)
        return NPOS;
    while ((baggage_type::eq(*(data()+count-1), c)) && (count > 1))
        --count;
    if ((count == 1) && (baggage_type::eq(*(data()), c)))
        return NPOS;
    else
        return count-1;
}

template <class charT>
basic_string<charT>
basic_string<charT>::substr (size_t pos,  size_t n) const _THROW_ALLOC_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if (length())
        return basic_string<charT> (data()+pos,
        (n > (length()-pos)) ? (length()-pos) : n);
    else
        return basic_string<charT>();
}

template <class charT>
int
basic_string<charT>::compare (const basic_string<charT>& str, size_t pos,
                              size_t n) const _THROW_OUTRANGE
{
    size_t slen   = (n > (length()-pos)) ? (length()-pos) : n;
    return compare_str (pos, str.data(), slen, str.length());
}

template <class charT>
int
basic_string<charT>::compare (const charT* s, size_t pos, size_t n) const
                              _THROW_LENGTH_OUTRANGE
{
    if (n == NPOS)
    {
        reference_class::throwlength();
    }
    return compare_str (pos, s, length()-pos, n);
}

template <class charT>
int
basic_string<charT>::compare (const charT* s, size_t pos) const _THROW_OUTRANGE
{
    return compare_str (pos, s, length()-pos, baggage_type::length(s));
}

template <class charT>
int
basic_string<charT>::compare (charT c, size_t pos, size_t rep) const
                              _THROW_LENGTH_OUTRANGE
{
    if (pos > length())
    {
        reference_class::throwrange();
    }
    if (rep == NPOS)
    {
        reference_class::throwlength();
    }
    if (rep)
    {
        size_t count = 0;
        while ((count < rep) && (count < (length()-pos)) &&
                baggage_type::eq (*(data()+pos+count), c))
            ++count;
        if ((count == rep) || (count == (length()-pos)))
            return (length()-pos-count);
        else
            return (*(data()+pos+count)-c);
    }
    else
    {
        return (length()-pos);
    }
}

template <class charT>
basic_string<charT>
operator+ (const basic_string<charT>& lhs, const basic_string<charT>& rhs)
           _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type;
    basic_string<charT> tmp(lhs.data(), lhs.length(), rhs.length());
    if (rhs.length())
        baggage_type::copy (tmp.point()+lhs.length(), rhs.data(), rhs.length());
    tmp.len() += rhs.length();
    return tmp;
}

template <class charT>
basic_string<charT>
operator+ (const charT* lhs, const basic_string<charT>& rhs) _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type; 
    size_t  slen = baggage_type::length(lhs);
    basic_string<charT> tmp(lhs, slen, rhs.length());
    if (rhs.length())
        baggage_type::copy (tmp.point()+slen, rhs.data(), rhs.length());
    tmp.len() += rhs.length();
    return tmp;
}

template <class charT>
basic_string<charT>
operator+ (charT lhs, const basic_string<charT>& rhs) _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type; 
    basic_string<charT> tmp(&lhs, 1, rhs.length());
    if (rhs.length())
        baggage_type::copy (tmp.point()+1, rhs.data(), rhs.length());
    tmp.len() += rhs.length();
    return tmp;
}

template <class charT>
basic_string<charT>
operator+ (const basic_string<charT>& lhs, const charT* rhs) _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type; 
    size_t  slen = baggage_type::length(rhs);
    basic_string<charT> tmp(lhs.data(), lhs.length(), slen);
    if (slen)
        baggage_type::copy (tmp.point()+lhs.length(), rhs, slen);
    tmp.len() += slen;
    return tmp;
}

template <class charT>
basic_string<charT>
operator+ (const basic_string<charT>& lhs, charT rhs) _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type; 
    basic_string<charT> tmp(lhs.data(), lhs.length(), 1);
    baggage_type::assign (*(tmp.point()+lhs.length()), rhs);
    ++tmp.len();
    return tmp;
}

template <class charT>
ostream&
operator<< (ostream& o, const basic_string<charT>& s) _THROW_NONE
{
    typedef  basic_string<charT>::baggage_type  baggage_type;
    for (size_t count = 0; count < s.length(); ++count)
        baggage_type::char_out (o, *(s.data()+count));
    return o;
}

template <class charT>
istream&
operator>> (istream& i, basic_string<charT>& s) _THROW_ALLOC_LENGTH
{
    typedef  basic_string<charT>::baggage_type  baggage_type; 
    s.remove();
    while (true)
    {
        charT  value;
        baggage_type::char_in (i, value);
        if (!i.operator void*())
            break;
        if (!baggage_type::is_del (value))
        {
            s.append(value);
            while (true)
            {
                baggage_type::char_in (i, value);
                if (!i.operator void*())
                    break;
                if (!baggage_type::is_del (value)) 
                {
                    s.append(value);
                }
                else
                    break;
            }
            break;
        }
    }
    return i;
}

typedef  basic_string<char>     cstring;
typedef  basic_string<char>     string;
// typedef  basic_string<wchar_t>  wstring;

struct hash<string>
{
  size_t operator()(const string& str) const
  {
    unsigned long h = 0;
    const string::char_type* s = str.data();
    for (size_t len = str.length(); len > 0; --len, ++s)
      h = 5*h + (unsigned long)(*s);
    return size_t(h);
  }
};

#if defined(__sgi) && !defined(__GNUC__) && (_MIPS_SIM != _MIPS_SIM_ABI32)
#pragma reset woff 1209
#endif


#endif /* __MBSTRING_H */
