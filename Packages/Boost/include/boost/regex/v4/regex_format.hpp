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
  *   FILE         regex_format.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Provides formatting output routines for search and replace
  *                operations.  Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */

#ifndef BOOST_REGEX_FORMAT_HPP
#define BOOST_REGEX_FORMAT_HPP


namespace boost{

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

//
// Forward declaration:
//
template <class RandomAccessIterator, class Allocator>
class match_results;

namespace re_detail{

// make_upper and make_lower should ideally be implemented in regex_traits
#if defined(_WIN32) && !defined(BOOST_REGEX_NO_W32)

//
// VC6 needs to link to user32.lib, as do all compilers that
// claim to be VC6/7 compatible:
//
#if defined(_MSC_VER) && !defined(__BORLANDC__)
#pragma comment(lib, "user32.lib")
#endif

inline wchar_t make_upper(wchar_t c)
{
   return LOWORD(::CharUpperW(reinterpret_cast<wchar_t*>(static_cast<unsigned short>(c))));
}

inline char make_upper(char c)
{
   return static_cast<char>(LOWORD(::CharUpperA(reinterpret_cast<char*>(static_cast<unsigned short>(c)))));
}

inline wchar_t make_lower(wchar_t c)
{
   return LOWORD(::CharLowerW(reinterpret_cast<wchar_t*>(static_cast<unsigned short>(c))));
}

inline char make_lower(char c)
{
   return static_cast<char>(LOWORD(::CharLowerA(reinterpret_cast<char*>(static_cast<unsigned short>(c)))));
}

#else

// TODO: make this traits class sensitive:
#ifndef BOOST_NO_WREGEX
inline wchar_t make_upper(wchar_t c)
{
   return (std::towupper)(c); 
}

inline wchar_t make_lower(wchar_t c)
{
   return (std::towlower)(c);
}
#endif
inline char make_upper(char c)
{
   return static_cast<char>((std::toupper)(c)); 
}

inline char make_lower(char c)
{
   return static_cast<char>((std::tolower)(c)); 
}

#endif //defined(_WIN32) && !defined(BOOST_REGEX_NO_W32)

typedef enum {
   case_nochange,
   case_oneupper,
   case_onelower,
   case_allupper,
   case_alllower
} case_flags_type;

// traits_type is unused, but provided to make it possible to use it for case conversion
template <class O, class charT, class traits_type>
void BOOST_REGEX_CALL output_char(O& out, charT c, traits_type& /*t*/, case_flags_type& f)
{
   switch (f) {
      case case_oneupper:
         f = case_nochange;
         // drop through
      case case_allupper:
         *out = make_upper(c);
         break;
      case case_onelower:
         f = case_nochange;
         // drop through
      case case_alllower:
         *out = make_lower(c);
         break;
      default:
         *out = c;
         break;
   }
}

template <class O, class I, class traits_type>
O BOOST_REGEX_CALL re_copy_out(O out, I first, I last, traits_type& t, case_flags_type& f)
{
   while(first != last)
   {
      if (f != case_nochange)
         output_char(out, *first, t, f);
      else
         *out = *first;

      ++out;
      ++first;
   }
   return out;
}


template <class charT, class traits_type>
void BOOST_REGEX_CALL re_skip_format(const charT*& fmt, const traits_type& traits_inst)
{
   // dwa 9/13/00 - suppress incorrect unused parameter warning for MSVC
   (void)traits_inst;

   typedef typename traits_type::size_type traits_size_type;
   typedef typename traits_type::uchar_type traits_uchar_type;
   typedef typename traits_type::string_type traits_string_type;

   unsigned int parens = 0;
   unsigned int c;
   while(*fmt)
   {
      c = traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt));
      if((c == traits_type::syntax_colon) && (parens == 0))
      {
         ++fmt;
         return;
      }
      else if(c == traits_type::syntax_close_bracket)
      {
         if(parens == 0)
         {
            ++fmt;
            return;
         }
         --parens;
      }
      else if(c == traits_type::syntax_open_bracket)
         ++parens;
      else if(c == traits_type::syntax_slash)
      {
         ++fmt;
         if(*fmt == 0)
            return;
      }
      ++fmt;
   }
}

#ifdef BOOST_NO_STD_OUTPUT_ITERATOR_ASSIGN

//
// ugly hack for buggy output iterators

template <class T>
inline void oi_assign(T* p, T v)
{
   ::boost::re_detail::pointer_destroy(p);
   pointer_construct(p, v);
}

#else

template <class T>
inline void oi_assign(T* p, T v)
{
   //
   // if you get a compile time error in here then you either
   // need to rewrite your output iterator to make it assignable
   // (as is required by the standard), or define
   // BOOST_NO_STD_OUTPUT_ITERATOR_ASSIGN to use the ugly hack above
   *p = v;
}

#endif


#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
//
// Ugly ugly hack,
// template don't merge if they contain switch statements so declare these
// templates in unnamed namespace (ie with internal linkage), each translation
// unit then gets its own local copy, it works seemlessly but bloats the app.
namespace{
#endif

//
// algorithm reg_format:
// takes the result of a match and a format string
// and merges them to produce a new string which
// is sent to an OutputIterator,
// _reg_format_aux does the actual work:
//
template <class OutputIterator, class Iterator, class Allocator, class charT, class traits_type>
OutputIterator BOOST_REGEX_CALL _reg_format_aux(OutputIterator out,
                          const match_results<Iterator, Allocator>& m,
                          const charT*& fmt,
                          match_flag_type flags, const traits_type& traits_inst,
                          case_flags_type& case_flags)
{
#ifdef __BORLANDC__
#pragma option push -w-8037
#endif
   const charT* fmt_end = fmt;
   while(*fmt_end) ++ fmt_end;

   typedef typename traits_type::size_type traits_size_type;
   typedef typename traits_type::uchar_type traits_uchar_type;
   typedef typename traits_type::string_type traits_string_type;

   while(*fmt)
   {
      switch(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)))
      {
      case traits_type::syntax_dollar:
         if(flags & format_sed)
         {
            // no perl style replacement,
            // $ is an ordinary character:
            goto default_opt;
         }
         ++fmt;
         if(*fmt == 0) // oops trailing $
         {
            --fmt;
            *out = *fmt;
            ++out;
            return out;
         }
         switch(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)))
         {
         case traits_type::syntax_start_buffer:
            oi_assign(&out, re_copy_out(out, Iterator(m[-1].first), Iterator(m[-1].second),
                        traits_inst, case_flags));
            ++fmt;
            continue;
         case traits_type::syntax_end_buffer:
            oi_assign(&out, re_copy_out(out, Iterator(m[-2].first), Iterator(m[-2].second),
                        traits_inst, case_flags));
            ++fmt;
            continue;
         case traits_type::syntax_digit:
         {
expand_sub:
            unsigned int index = traits_inst.toi(fmt, fmt_end, 10);
            if(index < m.size())
               oi_assign(&out, re_copy_out(out, Iterator(m[index].first), Iterator(m[index].second),
                        traits_inst, case_flags));
            continue;
         }
         }
         // anything else:
         if(*fmt == '&')
         {
            oi_assign(&out, re_copy_out(out, Iterator(m[0].first), Iterator(m[0].second),
                        traits_inst, case_flags));
            ++fmt;
         }
         else
         {
            // probably an error, treat as a literal '$'
            --fmt;
            *out = *fmt;
            ++out;
            ++fmt;
         }
         continue;
      case traits_type::syntax_slash:
      {
         // escape sequence:
         ++fmt;
         charT c(*fmt);
         if(*fmt == 0)
         {
            --fmt;
            *out = *fmt;
            ++out;
            ++fmt;
            return out;
         }
         switch(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)))
         {
         case traits_type::syntax_a:
            c = '\a';
            ++fmt;
            break;
         case traits_type::syntax_f:
            c = '\f';
            ++fmt;
            break;
         case traits_type::syntax_n:
            c = '\n';
            ++fmt;
            break;
         case traits_type::syntax_r:
            c = '\r';
            ++fmt;
            break;
         case traits_type::syntax_t:
            c = '\t';
            ++fmt;
            break;
         case traits_type::syntax_v:
            c = '\v';
            ++fmt;
            break;
         case traits_type::syntax_x:
            ++fmt;
            if(fmt == fmt_end)
            {
               *out = *--fmt;
               ++out;
               return out;
            }
            // maybe have \x{ddd}
            if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)) == traits_type::syntax_open_brace)
            {
               ++fmt;
               if(fmt == fmt_end)
               {
                  fmt -= 2;
                  *out = *fmt;
                  ++out;
                  ++fmt;
                  continue;
               }
               if(traits_inst.is_class(*fmt, traits_type::char_class_xdigit) == false)
               {
                  fmt -= 2;
                  *out = *fmt;
                  ++out;
                  ++fmt;
                  continue;
               }
               c = (charT)traits_inst.toi(fmt, fmt_end, -16);
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)) != traits_type::syntax_close_brace)
               {
                  while(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*fmt)) != traits_type::syntax_slash)
                     --fmt;
                  ++fmt;
                  *out = *fmt;
                  ++out;
                  ++fmt;
                  continue;
               }
               ++fmt;
               break;
            }
            else
            {
               if(traits_inst.is_class(*fmt, traits_type::char_class_xdigit) == false)
               {
                  --fmt;
                  *out = *fmt;
                  ++out;
                  ++fmt;
                  continue;
               }
               c = (charT)traits_inst.toi(fmt, fmt_end, -16);
            }
            break;
         case traits_type::syntax_c:
            ++fmt;
            if(fmt == fmt_end)
            {
               --fmt;
               *out = *fmt;
               ++out;
               return out;
            }
            if(((typename traits_type::uchar_type)(*fmt) < (typename traits_type::uchar_type)'@')
                  || ((typename traits_type::uchar_type)(*fmt) > (typename traits_type::uchar_type)127) )
            {
               --fmt;
               *out = *fmt;
               ++out;
               ++fmt;
               break;
            }
            c = (charT)((typename traits_type::uchar_type)(*fmt) - (typename traits_type::uchar_type)'@');
            ++fmt;
            break;
         case traits_type::syntax_e:
            c = (charT)27;
            ++fmt;
            break;
         case traits_type::syntax_digit:
            if(flags & format_sed)
               goto expand_sub;
            else
               c = (charT)traits_inst.toi(fmt, fmt_end, -8);
            break;

         case traits_type::syntax_u:
            ++fmt;
            if(flags & format_sed) break;
            case_flags = case_oneupper;
            continue;
         case traits_type::syntax_l:
            ++fmt;
            if(flags & format_sed) break;
            case_flags = case_onelower;
            continue;
         case traits_type::syntax_U:
            ++fmt;
            if(flags & format_sed) break;
            case_flags = case_allupper;
            continue;
         case traits_type::syntax_L:
            ++fmt;
            if(flags & format_sed) break;
            case_flags = case_alllower;
            continue;
         case traits_type::syntax_E:
            ++fmt;
            if(flags & format_sed) break;
            case_flags = case_nochange;
            continue;
         default:
            //c = *fmt;
            ++fmt;
         }
         *out = c;
         ++out;
         continue;
      }
      case traits_type::syntax_open_bracket:
         if(0 == (flags & format_all))
         {
            *out = *fmt;
            ++out;
            ++fmt;
            continue;
         }
         else
         {
            ++fmt;  // recurse
            oi_assign(&out, _reg_format_aux(out, m, fmt, flags, traits_inst, case_flags));
            continue;
         }
      case traits_type::syntax_close_bracket:
         if(0 == (flags & format_all))
         {
            *out = *fmt;
            ++out;
            ++fmt;
            continue;
         }
         else
         {
            ++fmt;  // return from recursion
            return out;
         }
      case traits_type::syntax_colon:
         if(flags & regex_constants::format_is_if)
         {
            ++fmt;
            return out;
         }
         *out = *fmt;
         ++out;
         ++fmt;
         continue;
      case traits_type::syntax_question:
      {
         if(0 == (flags & format_all))
         {
            *out = *fmt;
            ++out;
            ++fmt;
            continue;
         }
         else
         {
            ++fmt;
            if(*fmt == 0)
            {
               --fmt;
               *out = *fmt;
               ++out;
               ++fmt;
               return out;
            }
            unsigned int id = traits_inst.toi(fmt, fmt_end, 10);
            if(m[id].matched)
            {
               oi_assign(&out, _reg_format_aux(out, m, fmt, flags | regex_constants::format_is_if, traits_inst, case_flags));
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*(fmt-1))) == traits_type::syntax_colon)
                  re_skip_format(fmt, traits_inst);
            }
            else
            {
               re_skip_format(fmt, traits_inst);
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*(fmt-1))) == traits_type::syntax_colon)
                  oi_assign(&out, _reg_format_aux(out, m, fmt, flags | regex_constants::format_is_if, traits_inst, case_flags));
            }
            return out;
         }
      }
      default:
default_opt:
         if((flags & format_sed) && (*fmt == '&'))
         {
            oi_assign(&out, re_copy_out(out, Iterator(m[0].first), Iterator(m[0].second),
                     traits_inst, case_flags));
            ++fmt;
            continue;
         }

         output_char(out, *fmt, traits_inst, case_flags);
         ++out;
         ++fmt;
      }
   }

   return out;
#ifdef __BORLANDC__
#pragma option pop
#endif
}

#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
} // namespace
#endif

template <class S>
class string_out_iterator
{
   S* out;
public:
   typedef typename S::difference_type difference_type;
   typedef typename S::value_type value_type;
   typedef typename S::pointer pointer;
   typedef typename S::reference reference;
   typedef          std::output_iterator_tag iterator_category;

   string_out_iterator(S& s) : out(&s) {}
   string_out_iterator& operator++() { return *this; }
   string_out_iterator& operator++(int) { return *this; }
   string_out_iterator& operator*() { return *this; }
   string_out_iterator& operator=(typename S::value_type v)
   {
      out->append(1, v);
      return *this;
   }
};

template <class OutputIterator, class Iterator, class charT, class Allocator, class traits_type>
class merge_out_predicate
{
   OutputIterator* out;
   Iterator* last;
   const charT* fmt;
   match_flag_type flags;
   const traits_type* pt;
   // rebind allocator to correct type:
   typedef typename detail::rebind_allocator<sub_match<Iterator>, Allocator>::type alloc_type;

public:
   merge_out_predicate(OutputIterator& o, Iterator& pi, const charT* f, match_flag_type format_flags, const traits_type& p)
      : out(&o), last(&pi), fmt(f), flags(format_flags), pt(&p){}

   ~merge_out_predicate() {}
   bool BOOST_REGEX_CALL operator()(const boost::match_results<Iterator, alloc_type>& m)
   {
      const charT* f = fmt;
      case_flags_type cf = case_nochange;
      if(0 == (flags & format_no_copy))
      {
         oi_assign(out, re_copy_out(
                     *out, 
                     Iterator(m[-1].first), 
                     Iterator(m[-1].second),
                     *pt, 
                     cf));
      }
      oi_assign(out, _reg_format_aux(*out, m, f, flags, *pt, cf));
      *last = m[-2].first;
      return flags & format_first_only ? false : true;
   }
};

} // namespace re_detail

template <class OutputIterator, class Iterator, class Allocator, class charT>
OutputIterator regex_format(OutputIterator out,
                          const match_results<Iterator, Allocator>& m,
                          const charT* fmt,
                          match_flag_type flags = format_all
                         )
{
   regex_traits<charT> t;

   re_detail::case_flags_type cf = re_detail::case_nochange;
   return re_detail::_reg_format_aux(out, m, fmt, flags, t, cf);
}

template <class OutputIterator, class Iterator, class Allocator, class charT>
OutputIterator regex_format(OutputIterator out,
                          const match_results<Iterator, Allocator>& m,
                          const std::basic_string<charT>& fmt,
                          match_flag_type flags = format_all
                         )
{
   regex_traits<charT> t;
   const charT* start = fmt.c_str();

   re_detail::case_flags_type cf = re_detail::case_nochange;
   return re_detail::_reg_format_aux(out, m, start, flags, t, cf);
}

template <class Iterator, class Allocator, class charT>
std::basic_string<charT> regex_format(const match_results<Iterator, Allocator>& m,
                                      const charT* fmt,
                                      match_flag_type flags = format_all)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_format(i, m, fmt, flags);
   return result;
}

template <class Iterator, class Allocator, class charT>
std::basic_string<charT> regex_format(const match_results<Iterator, Allocator>& m,
                                      const std::basic_string<charT>& fmt,
                                      match_flag_type flags = format_all)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_format(i, m, fmt.c_str(), flags);
   return result;
}

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost

#endif  // BOOST_REGEX_FORMAT_HPP






