// -*- C++ -*-
//  Boost general library 'format'   ---------------------------
//  See http://www.boost.org for updates, documentation, and revision history.

//  (C) Samuel Krempp 2001
//                  krempp@crans.ens-cachan.fr
//  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

// ideas taken from Rüdiger Loos's format class
// and Karl Nelson's ofstream (also took its parsing code as basis for printf parsing)

// ------------------------------------------------------------------------------
// format_class.hpp :  class interface
// ------------------------------------------------------------------------------

#ifndef BOOST_FORMAT_CLASS_HPP
#define BOOST_FORMAT_CLASS_HPP


#include <vector>
#include <string>

#include <boost/format/format_fwd.hpp>
#include <boost/format/internals_fwd.hpp>
#include <boost/format/internals.hpp>
#include <boost/format/outsstream.hpp>

namespace boost {

    template<class Ch, class Tr>
    class basic_format 
    {
    public:
        typedef Ch  CharT;   // borland fails if we use Ch and Tr directly
        typedef Tr  Traits;  // in the body of the operator% template.
        typedef std::basic_string<Ch, Tr>                string_t;
        typedef typename std::basic_string<Ch, Tr>::size_type size_type; 

        typedef io::basic_outsstream<Ch, Tr>             internal_stream_t;

        explicit basic_format(const Ch* str=NULL);
        explicit basic_format(const string_t& s);
#ifndef BOOST_NO_STD_LOCALE
        explicit basic_format(const Ch* str, const std::locale & loc);
        explicit basic_format(const string_t& s, const std::locale & loc);
#endif
        basic_format(const basic_format& x);
        basic_format& operator= (const basic_format& x);
        void swap(basic_format& x);

        basic_format& clear(); // empty all converted string buffers (except bound items)
        basic_format& clear_binds(); // unbind all bound items, and call clear()
        basic_format& parse(const string_t&); // resets buffers and parse a new format string

        // pass arguments through those operators :
        template<class T>  
        basic_format&   operator%(const T& x)
            { return io::detail::feed<CharT, Traits, const T&>(*this,x); }

#ifndef BOOST_NO_OVERLOAD_FOR_NON_CONST
        template<class T>  basic_format&   operator%(T& x) 
            { return io::detail::feed<CharT, Traits, T&>(*this,x); }
#endif
        // modifying a format object
        template<class T>
        basic_format&  bind_arg(int argN, const T& val) 
            { return io::detail::bind_arg_body(*this, argN, val); }
        basic_format&  clear_bind(int argN);
        template<class T> 
        basic_format&  modify_item(int itemN, T manipulator) 
            { return io::detail::modify_item_body<Ch,Tr,T> (*this, itemN, manipulator);}

        // Choosing which errors will throw exceptions :
        unsigned char exceptions() const;
        unsigned char exceptions(unsigned char newexcept);

        // final output
        size_type size() const;            // sum of the current string pieces sizes
        string_t str() const;
        friend BOOST_IO_STD basic_ostream<Ch, Tr>& 
#if BOOST_WORKAROUND( BOOST_MSVC, <= 1300) 
        operator<< (BOOST_IO_STD basic_ostream<Ch, Tr>& , const basic_format& ); 
#else
        operator<< <Ch, Tr> (BOOST_IO_STD basic_ostream<Ch, Tr>&, const basic_format&); 
#endif
                      
#if !defined( BOOST_NO_MEMBER_TEMPLATE_FRIENDS )  \
    && !BOOST_WORKAROUND(__BORLANDC__, <= 0x570) \
    && !BOOST_WORKAROUND( _CRAYC, != 0)
        // use friend templates and private members only if supported
        template<class Ch2, class Tr2, class T>  friend basic_format<Ch2, Tr2>&  
        io::detail::feed(basic_format<Ch2,Tr2>&, T);
        template<class Ch2, class Tr2, class T>  friend   
        void io::detail::distribute(basic_format<Ch2,Tr2>&, T);
        template<class Ch2, class Tr2, class T>  friend
        basic_format<Ch2, Tr2>& io::detail::modify_item_body(basic_format<Ch2, Tr2>&, int, T);
        template<class Ch2, class Tr2, class T> friend
        basic_format<Ch2, Tr2>&  io::detail::bind_arg_body(basic_format<Ch2, Tr2>&, int, const T&);

    private:
#endif
        typedef BOOST_IO_STD basic_ostream<Ch, Tr>       stream_t;
        typedef io::detail::stream_format_state<Ch, Tr>  stream_format_state;
        typedef io::detail::format_item<Ch, Tr>          format_item_t;
        // flag bits, used for style_
        enum style_values  { ordered = 1, // set only if all directives are  positional
                             special_needs = 4 };     

        void make_or_reuse_data(std::size_t nbitems);// used for (re-)initialisation

        std::vector<format_item_t>  items_; // each '%..' directive leads to a format_item
        std::vector<bool> bound_; // stores which arguments were bound. size() == 0 || num_args
        int               style_; // style of format-string :  positional or not, etc
        int             cur_arg_; // keep track of wich argument is current
        int            num_args_; // number of expected arguments
        mutable bool     dumped_; // true only after call to str() or <<
        string_t         prefix_; // piece of string to insert before first item
        internal_stream_t   oss_; // the internal stream.
        unsigned char exceptions_;
    }; // class basic_format

} // namespace boost


#endif // BOOST_FORMAT_CLASS_HPP
