// -*- C++ -*-
//  Boost  format   ----------------------------------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2003
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ------------------------------------------------------------------------------
// implementation included by outsstream.hpp.
// ------------------------------------------------------------------------------

namespace boost {
namespace io {

#if !defined(BOOST_NO_STRINGSTREAM) && !defined(BOOST_FORMAT_IGNORE_STRINGSTREAM)

template<class Ch, class Tr> inline
void steal_basic_stringbuf<Ch, Tr> :: clear_buffer() { 
    const Ch * p = pptr();
    const Ch * b = pbase();
    if(p != NULL && p != b) {
      typedef typename Tr::pos_type pos_type;
      pos_type pos = buff_t::seekpos(0, std::ios_base::out); 
      BOOST_ASSERT( pos != pos_type(std::streamoff(-1)) ); 
    }
}



#else // BOOST_NO_STRINGSTREAM


template <class Tr> inline
basic_outsstream<char,Tr> ::basic_outsstream() : 
  buff_t(),  
  stream_t(this) 
{ 
  stream_t::init(this);  // the strem construction isnt enough with gcc-2.95
}  
  

template <class Tr>
std::basic_string<char, Tr> 
basic_outsstream<char, Tr> ::str() {   // ! physically copy chars :
    string_type s(buff_t::str(), pcount());
    freeze(false);
    return s;
}

template<class Tr > inline
void
basic_outsstream<char, Tr>:: clear_buffer() { 
    freeze(false);
    const Ch * p = cur();
    const Ch * b = begin();
    if(p != NULL && p != b) {
      buff_t::seekpos(0, std::ios_base::out); 
    }
    freeze(false);
}

#endif

} //namespace io
} //namespace boost
