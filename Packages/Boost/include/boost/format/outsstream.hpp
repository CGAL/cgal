// -*- C++ -*-
//  Boost  format   ----------------------------------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2003
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ------------------------------------------------------------------------------
//   outsstream<Ch, Tr>   is  a class extending stringstream by adding :
//     clear_buffer() method, and
//     access to the [pbase(), pptr()[ part of the 'output sequence' 
//     (see §27.5.1 of the C++ standard)
//     if sstream is not available, it uses strstream and is slightly different.
// ------------------------------------------------------------------------------

#ifndef BOOST_FORMAT_OUTSSTREAM_H
#define BOOST_FORMAT_OUTSSTREAM_H

#include <boost/assert.hpp>
#include <boost/utility/base_from_member.hpp>
#include <boost/format/detail/config_macros.hpp>

#if !defined(BOOST_NO_STRINGSTREAM) && !defined(BOOST_FORMAT_IGNORE_STRINGSTREAM)
#include <sstream>
#else
#include <strstream>
#include <string>
#endif // BOOST_NO_STRING_STREAM


namespace boost {
namespace io {

#if !defined(BOOST_NO_STRINGSTREAM) && !defined(BOOST_FORMAT_IGNORE_STRINGSTREAM)


//---- The stringstream way ---------------------------------------------------


// ---   class steal_basic_stringbuf : steals  pbase(), pptr() & co -----------
    template<class Ch, class Tr = BOOST_IO_STD char_traits<Ch> >
    class steal_basic_stringbuf : public std::basic_stringbuf<Ch, Tr>
    {
        typedef std::basic_stringbuf<Ch, Tr> buff_t;
    public:
        typedef std::basic_string<Ch,Tr>     string_type;

        // get [pbase, pptr[   from stringbuf::str(),  which returns [pbase, epptr[  :
        string_type cur_str() const { return string_type(this->str(), 0, pcount()); } 

        // publicize those functions (protected in stringbuf) :
        std::streamsize pcount() const { return pptr() - pbase(); }
        const Ch * pbase() const { return buff_t::pbase(); } 
        const Ch * pptr()  const { return buff_t::pptr(); } 
        const Ch * epptr() const { return buff_t::epptr(); }
        // added convenience function :
        void clear_buffer();
    };


// ---   class basic_outsstream -----------------------------------------------
    template<class Ch, class Tr = BOOST_IO_STD char_traits<Ch> >
    class basic_outsstream : boost::base_from_member<steal_basic_stringbuf<Ch, Tr> >, 
                             public BOOST_IO_STD basic_ostream<Ch, Tr>
    // use stringbuf with its stolen protected members, 
    // publicly derived from basic_ostream to make this class a stream.
    {
    public:
        typedef std::basic_string<Ch,Tr>     string_type;
        basic_outsstream() : pbase_type(),
                             std::basic_ostream<Ch,Tr>( & sb() ) {}
        // buffer access via strings
        string_type str()     const { return sb().str(); }     // [begin, end[
        string_type cur_str() const { return sb().cur_str(); } // [begin, cur[
            
        // copy-less access (note the pointers can be invalidated when modifying the stream)
        std::streamsize pcount() const { return sb().pcount(); }
        const Ch * begin() const { return sb().pbase(); } 
        const Ch * cur()   const { return sb().pptr(); } 
        const Ch * end()   const { return sb().epptr(); }

        void clear_buffer() { sb().clear_buffer(); }
    private:
        typedef boost::base_from_member<steal_basic_stringbuf<Ch, Tr> > pbase_type;
        steal_basic_stringbuf<Ch, Tr>      &  sb()       { return pbase_type::member; }
        steal_basic_stringbuf<Ch, Tr> const&  sb() const { return pbase_type::member; }
    };
#else // BOOST_NO_STRINGSTREAM



//---- The strstream way ------------------------------------------------------

    template <class Ch, 
#if !( BOOST_WORKAROUND(__GNUC__, <3) &&  defined(__STL_CONFIG_H) )
        class Tr = BOOST_IO_STD char_traits<Ch> > 
#else
    class Tr = std::string_char_traits<Ch> > 
#endif
    class basic_outsstream; // we define only the <char> specialisaton


// ---   class basic_outsstream -----------------------------------------------
    template<class Tr>
    class basic_outsstream<char, Tr> : private BOOST_IO_STD  strstreambuf, 
                                       public std::basic_ostream<char, Tr>
    {
    public:
        typedef char Ch;
        typedef std::basic_string<Ch,Tr>      string_type;
        typedef BOOST_IO_STD strstreambuf     buff_t;
        typedef std::basic_ostream<char, Tr>  stream_t;
    public:
        basic_outsstream();

        // ! physically copy chars :
        string_type str(); 
        string_type cur_str() const  { return string_type(begin(), cur());  } 

        int freeze() const { return buff_t::freeze(); }
        void freeze(int x) { buff_t::freeze(x); }

        // copy-less access (be careful, the pointer are invalidated when modifying the stream) :
        std::streamsize pcount() const { return cur()-begin(); }
        const Ch * begin() const { return buff_t::pbase(); } 
        const Ch * cur()   const { return buff_t::pptr(); } 
        const Ch * end()   const { return buff_t::epptr(); }

        void clear_buffer();
    };
  
#endif // BOOST_NO_STRINGSTREAM


    typedef basic_outsstream<char> outsstream;

} // namespace boost::io
} // namespace boost


#include <boost/format/outsstream_impl.hpp> // implementation

#endif // BOOST_FORMAT_OUTSSTREAM_H include guard
