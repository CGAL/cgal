/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         basic_regex.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Declares template class basic_regex (note that member function
  *                bodies are in regex_compile.hpp).
  */

#ifndef BOOST_REGEX_V4_BASIC_REGEX_HPP
#define BOOST_REGEX_V4_BASIC_REGEX_HPP

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost{
//
// class reg_expression
// represents the compiled
// regular expression:
//

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable : 4251 4231 4660)
#endif

#ifdef BOOST_REGEX_NO_FWD
template <class charT, class traits = regex_traits<charT>, class Allocator = BOOST_DEFAULT_ALLOCATOR(charT) >
#else
template <class charT, class traits, class Allocator >
#endif
class reg_expression : public regbase
{
public:
   // typedefs:
   typedef typename traits::size_type            traits_size_type;
   typedef typename traits::uchar_type           traits_uchar_type;
   typedef typename traits::string_type          traits_string_type;
   typedef charT                                 char_type;
   typedef traits                                traits_type;

   typedef charT                                 value_type;
   typedef charT&                                reference;
   typedef const charT&                          const_reference;
   typedef const charT*                          const_iterator;
   typedef const_iterator                        iterator;
   typedef typename Allocator::difference_type   difference_type;
   typedef typename Allocator::size_type         size_type;   
   typedef Allocator                             allocator_type;
   typedef Allocator                             alloc_type;
   typedef regex_constants::syntax_option_type   flag_type;
   // locale_type
   // placeholder for actual locale type used by the
   // traits class to localise *this.
   typedef typename traits::locale_type          locale_type;
   
public:
   explicit reg_expression(const Allocator& a = Allocator());
   explicit reg_expression(const charT* p, flag_type f = regex_constants::normal, const Allocator& a = Allocator());
   reg_expression(const charT* p1, const charT* p2, flag_type f = regex_constants::normal, const Allocator& a = Allocator());
   reg_expression(const charT* p, size_type len, flag_type f, const Allocator& a = Allocator());
   reg_expression(const reg_expression&);
   ~reg_expression();
   reg_expression& BOOST_REGEX_CALL operator=(const reg_expression&);
   reg_expression& BOOST_REGEX_CALL operator=(const charT* ptr)
   {
      set_expression(ptr, regex_constants::normal | regex_constants::use_except);
      return *this;
   }

   //
   // assign:
   reg_expression& assign(const reg_expression& that)
   { return *this = that; }
   reg_expression& assign(const charT* ptr, flag_type f = regex_constants::normal)
   {
      set_expression(ptr, f | regex_constants::use_except);
      return *this;
   }
   reg_expression& assign(const charT* ptr, size_type len, flag_type f)
   {
      std::basic_string<charT> s(ptr, len);
      set_expression(s.c_str(), f | regex_constants::use_except);
      return *this;
   }

   reg_expression& assign(const charT* arg_first,
                          const charT* arg_last,
                          flag_type f = regex_constants::normal)
   {
      set_expression(arg_first, arg_last, f | regex_constants::use_except);
      return *this;
   }
#if !defined(BOOST_NO_MEMBER_TEMPLATES) && !defined(__IBMCPP__)

   template <class ST, class SA>
   unsigned int BOOST_REGEX_CALL set_expression(const std::basic_string<charT, ST, SA>& p, flag_type f = regex_constants::normal)
   { return set_expression(p.data(), p.data() + p.size(), f); }

   template <class ST, class SA>
   explicit reg_expression(const std::basic_string<charT, ST, SA>& p, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
    : data(a), pkmp(0), error_code_(REG_EMPTY), _expression(0) { set_expression(p, f | regex_constants::use_except); }

   template <class InputIterator>
   reg_expression(InputIterator arg_first, InputIterator arg_last, flag_type f = regex_constants::normal, const Allocator& al = Allocator())
    : data(al), pkmp(0), error_code_(REG_EMPTY), _expression(0)
   {
      std::basic_string<charT> a(arg_first, arg_last);
      set_expression(a.data(), a.data() + a.size(), f | regex_constants::use_except);
   }

   template <class ST, class SA>
   reg_expression& BOOST_REGEX_CALL operator=(const std::basic_string<charT, ST, SA>& p)
   {
      set_expression(p.c_str(), p.c_str() + p.size(), regex_constants::normal | regex_constants::use_except);
      return *this;
   }

   template <class string_traits, class A>
   reg_expression& BOOST_REGEX_CALL assign(
       const std::basic_string<charT, string_traits, A>& s,
       flag_type f = regex_constants::normal)
   {
      set_expression(s.c_str(), s.c_str() + s.size(), f | regex_constants::use_except);
      return *this;
   }

   template <class InputIterator>
   reg_expression& BOOST_REGEX_CALL assign(InputIterator arg_first,
                          InputIterator arg_last,
                          flag_type f = regex_constants::normal)
   {
      std::basic_string<charT> a(arg_first, arg_last);
      set_expression(a.data(), a.data() + a.size(), f | regex_constants::use_except);
      return *this;
   }
#else
   unsigned int BOOST_REGEX_CALL set_expression(const std::basic_string<charT>& p, flag_type f = regex_constants::normal)
   { return set_expression(p.data(), p.data() + p.size(), f | regex_constants::use_except); }

   reg_expression(const std::basic_string<charT>& p, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
    : data(a), pkmp(0) { set_expression(p, f | regex_constants::use_except); }

   reg_expression& BOOST_REGEX_CALL operator=(const std::basic_string<charT>& p)
   {
      set_expression(p.c_str(), p.c_str() + p.size(), regex_constants::normal | regex_constants::use_except);
      return *this;
   }

   reg_expression& BOOST_REGEX_CALL assign(
       const std::basic_string<charT>& s,
       flag_type f = regex_constants::normal)
   {
      set_expression(s.c_str(), s.c_str() + s.size(), f | regex_constants::use_except);
      return *this;
   }

#endif


   //
   // allocator access:
   Allocator BOOST_REGEX_CALL get_allocator()const;
   //
   // locale:
   locale_type BOOST_REGEX_CALL imbue(locale_type l){ return traits_inst.imbue(l); }
   locale_type BOOST_REGEX_CALL getloc()const{ return traits_inst.getloc(); }
   //
   // getflags:
   // retained for backwards compatibility only, the base class has "flags"
   // member which is now the prefered name:
   flag_type BOOST_REGEX_CALL getflags()const
   { return this->flags(); }
   //
   // str:
   std::basic_string<charT> BOOST_REGEX_CALL str()const
   {
      std::basic_string<charT> result;
      if(this->error_code() == 0)
         result = std::basic_string<charT>(_expression, _expression_len);
      return result;
   }
   //
   // begin, end:
   const_iterator BOOST_REGEX_CALL begin()const
   { return (this->error_code() ? 0 : _expression); }
   const_iterator BOOST_REGEX_CALL end()const
   { return (this->error_code() ? 0 : _expression + _expression_len); }
   //
   // swap:
   void BOOST_REGEX_CALL swap(reg_expression&)throw();
   //
   // size:
   size_type BOOST_REGEX_CALL size()const
   { return (this->error_code() ? 0 : _expression_len); }
   //
   // max_size:
   size_type BOOST_REGEX_CALL max_size()const
   { return UINT_MAX; }
   //
   // empty:
   bool BOOST_REGEX_CALL empty()const
   { return 0 != this->error_code(); }

   unsigned BOOST_REGEX_CALL mark_count()const { return (this->error_code() ? 0 : marks); }
   int BOOST_REGEX_CALL compare(const reg_expression&) const;
   bool BOOST_REGEX_CALL operator==(const reg_expression& e)const
   { return compare(e) == 0; }
   bool operator != (const basic_regex<charT, traits, Allocator>& e)
   { return compare(e) != 0; }
   bool BOOST_REGEX_CALL operator<(const reg_expression& e)const
   { return compare(e) < 0; }
   bool BOOST_REGEX_CALL operator>(const reg_expression& e)const
   { return compare(e) > 0; }
   bool BOOST_REGEX_CALL operator<=(const reg_expression& e)const
   { return compare(e) <= 0; }
   bool BOOST_REGEX_CALL operator>=(const reg_expression& e)const
   { return compare(e) >= 0; }

   //
   // The following are deprecated as public interfaces
   // but are available for compatibility with earlier versions.
   allocator_type BOOST_REGEX_CALL allocator()const;
   const charT* BOOST_REGEX_CALL expression()const { return (this->error_code() ? 0 : _expression); }
   unsigned int BOOST_REGEX_CALL set_expression(const charT* p, const charT* end, flag_type f = regex_constants::normal);
   unsigned int BOOST_REGEX_CALL set_expression(const charT* p, flag_type f = regex_constants::normal) { return set_expression(p, p + traits_type::length(p), f); }
   //
   // this should be private but template friends don't work:
   const traits_type& get_traits()const { return traits_inst; }
   unsigned int BOOST_REGEX_CALL error_code()const
   {
      return error_code_;
   }

private:
   traits_type                          traits_inst;          // traits class in use
   re_detail::raw_storage<Allocator>    data;                 // our state machine
   unsigned                             _restart_type;        // search method to use
   unsigned                             marks;                // number of marked sub-expressions
   int                                  repeats;              // number of repeats
   unsigned char*                       startmap;             // characters that can match the first state(s) in the machine
   std::size_t                          _expression_len;      // length of the expression
   std::size_t                          _leading_len;         // length of any leading literal 
   const charT*                         _leading_string;      // leading literal string
   std::size_t                          _leading_string_len;  // and it's length
   re_detail::kmp_info<charT>*          pkmp;                 // pointer to Knuth Morris Pratt state machine when available
   unsigned                             error_code_;          // our current status
   charT*                               _expression;          // the expression we just compiled if any

   void BOOST_REGEX_CALL compile_maps();
   void BOOST_REGEX_CALL compile_map(re_detail::re_syntax_base* node, unsigned char* _map, unsigned int* pnull, unsigned char mask, re_detail::re_syntax_base* terminal = 0)const;
   bool BOOST_REGEX_CALL probe_start(re_detail::re_syntax_base* node, charT c, re_detail::re_syntax_base* terminal)const;
   bool BOOST_REGEX_CALL probe_start_null(re_detail::re_syntax_base* node, re_detail::re_syntax_base* terminal)const;
   void BOOST_REGEX_CALL fixup_apply(re_detail::re_syntax_base* b, unsigned cbraces);
   void BOOST_REGEX_CALL move_offsets(re_detail::re_syntax_base* j, unsigned size);
   re_detail::re_syntax_base* BOOST_REGEX_CALL compile_set(const charT*& first, const charT* last);
   re_detail::re_syntax_base* BOOST_REGEX_CALL compile_set_aux(re_detail::jstack<traits_string_type, Allocator>& singles, re_detail::jstack<traits_string_type, Allocator>& ranges, re_detail::jstack<boost::uint_fast32_t, Allocator>& classes, re_detail::jstack<traits_string_type, Allocator>& equivalents, bool isnot, const re_detail::_narrow_type&);
   re_detail::re_syntax_base* BOOST_REGEX_CALL compile_set_aux(re_detail::jstack<traits_string_type, Allocator>& singles, re_detail::jstack<traits_string_type, Allocator>& ranges, re_detail::jstack<boost::uint_fast32_t, Allocator>& classes, re_detail::jstack<traits_string_type, Allocator>& equivalents, bool isnot, const re_detail::_wide_type&);
   re_detail::re_syntax_base* BOOST_REGEX_CALL compile_set_simple(re_detail::re_syntax_base* dat, unsigned long cls, bool isnot = false);
   unsigned int BOOST_REGEX_CALL parse_inner_set(const charT*& first, const charT* last);

   re_detail::re_syntax_base* BOOST_REGEX_CALL add_simple(re_detail::re_syntax_base* dat, re_detail::syntax_element_type type, unsigned int size = sizeof(re_detail::re_syntax_base));
   re_detail::re_syntax_base* BOOST_REGEX_CALL add_literal(re_detail::re_syntax_base* dat, charT c);
   charT BOOST_REGEX_CALL parse_escape(const charT*& first, const charT* last);
   void BOOST_REGEX_CALL parse_range(const charT*& first, const charT* last, unsigned& min, unsigned& max);
   bool BOOST_REGEX_CALL skip_space(const charT*& first, const charT* last);
   unsigned int BOOST_REGEX_CALL probe_restart(re_detail::re_syntax_base* dat);
   unsigned int BOOST_REGEX_CALL fixup_leading_rep(re_detail::re_syntax_base* dat, re_detail::re_syntax_base* end);
   void BOOST_REGEX_CALL fail(unsigned int err);

protected:
   static int BOOST_REGEX_CALL repeat_count(const reg_expression& e)
   { return e.repeats; }
   static unsigned int BOOST_REGEX_CALL restart_type(const reg_expression& e)
   { return e._restart_type; }
   static const re_detail::re_syntax_base* BOOST_REGEX_CALL first(const reg_expression& e)
   { return (const re_detail::re_syntax_base*)e.data.data(); }
   static const unsigned char* BOOST_REGEX_CALL get_map(const reg_expression& e)
   { return e.startmap; }
   static std::size_t BOOST_REGEX_CALL leading_length(const reg_expression& e)
   { return e._leading_len; }
   static const re_detail::kmp_info<charT>* get_kmp(const reg_expression& e)
   { return e.pkmp; }
   static bool BOOST_REGEX_CALL can_start(charT c, const unsigned char* _map, unsigned char mask, const re_detail::_wide_type&);
   static bool BOOST_REGEX_CALL can_start(charT c, const unsigned char* _map, unsigned char mask, const re_detail::_narrow_type&);
};

template <class charT, class traits, class Allocator >
void swap(reg_expression<charT, traits, Allocator>& a, reg_expression<charT, traits, Allocator>& b)
{
   a.swap(b);
}

#ifndef BOOST_NO_STD_LOCALE
template <class charT, class traits, class traits2, class Allocator>
std::basic_ostream<charT, traits>& 
   operator << (std::basic_ostream<charT, traits>& os, 
                const reg_expression<charT, traits2, Allocator>& e)
{
   return (os << e.str());
}
#else
template <class traits, class Allocator>
std::ostream& operator << (std::ostream& os, const reg_expression<char, traits, Allocator>& e)
{
   return (os << e.str());
}
#endif

//
// We want to rename reg_expression basic_regex but maintain 
// backwards compatibility, so class basic_regex is just a thin
// wrapper around reg_expression:
//
#ifdef BOOST_REGEX_NO_FWD
template <class charT, class traits = regex_traits<charT>, class Allocator = BOOST_DEFAULT_ALLOCATOR(charT) >
#else
template <class charT, class traits, class Allocator >
#endif
class basic_regex : public reg_expression<charT, traits, Allocator>
{
public:
   typedef typename reg_expression<charT, traits, Allocator>::flag_type flag_type;
   typedef typename reg_expression<charT, traits, Allocator>::size_type size_type;
   explicit basic_regex(const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(a){}
   explicit basic_regex(const charT* p, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(p,f,a){}
   basic_regex(const charT* p1, const charT* p2, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(p1,p2,f,a){}
   basic_regex(const charT* p, size_type len, flag_type f, const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(p,len,f,a){}
   basic_regex(const basic_regex& that)
      : reg_expression<charT, traits, Allocator>(that){}
   ~basic_regex(){}
   basic_regex& BOOST_REGEX_CALL operator=(const basic_regex& that)
   {
      this->assign(that);
      return *this;
   }
   basic_regex& BOOST_REGEX_CALL operator=(const charT* ptr)
   {
      this->assign(ptr);
      return *this;
   }
#if !defined(BOOST_NO_MEMBER_TEMPLATES) && !(defined(__IBMCPP__) && (__IBMCPP__ <= 502))
   template <class ST, class SA>
   explicit basic_regex(const std::basic_string<charT, ST, SA>& p, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(p,f,a){}

   template <class I>
   basic_regex(I arg_first, I arg_last, flag_type f = regex_constants::normal, const Allocator& al = Allocator())
      : reg_expression<charT, traits, Allocator>(arg_first, arg_last, f, al){}

   template <class ST, class SA>
   basic_regex& BOOST_REGEX_CALL operator=(const std::basic_string<charT, ST, SA>& p)
   {
      this->assign(p);
      return *this;
   }
#else
   basic_regex(const std::basic_string<charT>& p, flag_type f = regex_constants::normal, const Allocator& a = Allocator())
      : reg_expression<charT, traits, Allocator>(p,f,a){}

   basic_regex& BOOST_REGEX_CALL operator=(const std::basic_string<charT>& p)
   {
      this->assign(p);
      return *this;
   }
#endif
};

#ifdef BOOST_MSVC
#pragma warning (pop)
#endif

} // namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif

