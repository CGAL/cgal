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

enum format_flags_t{
   format_all = 0,                      // enable all extentions to sytax
   format_sed = match_max << 1,         // sed style replacement.
   format_perl = format_sed << 1,       // perl style replacement.
   format_no_copy = format_perl << 1,   // don't copy non-matching segments.
   format_first_only = format_no_copy << 1,   // Only replace first occurance.
   format_is_if = format_first_only << 1   // internal use only.
};

namespace re_detail{

#ifdef __BORLANDC__
   #pragma option push -a8 -b -Vx -Ve -pc  -w-8037
#endif

template <class O, class I>
O BOOST_REGEX_CALL re_copy_out(O out, I first, I last)
{
   while(first != last)
   {
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
                          unsigned flags, const traits_type& traits_inst)
{
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
            oi_assign(&out, re_copy_out(out, Iterator(m[-1].first), Iterator(m[-1].second)));
            ++fmt;
            continue;
         case traits_type::syntax_end_buffer:
            oi_assign(&out, re_copy_out(out, Iterator(m[-2].first), Iterator(m[-2].second)));
            ++fmt;
            continue;
         case traits_type::syntax_digit:
         {
expand_sub:
            unsigned int index = traits_inst.toi(fmt, fmt_end, 10);
            oi_assign(&out, re_copy_out(out, Iterator(m[index].first), Iterator(m[index].second)));
            continue;
         }
         }
         // anything else:
         if(*fmt == '&')
         {
            oi_assign(&out, re_copy_out(out, Iterator(m[0].first), Iterator(m[0].second)));
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
         default:
            //c = *fmt;
            ++fmt;
         }
         *out = c;
         ++out;
         continue;
      }
      case traits_type::syntax_open_bracket:
         if(flags & (format_sed|format_perl))
         {
            *out = *fmt;
            ++out;
            ++fmt;
            continue;
         }
         else
         {
            ++fmt;  // recurse
            oi_assign(&out, _reg_format_aux(out, m, fmt, flags, traits_inst));
            continue;
         }
      case traits_type::syntax_close_bracket:
         if(flags & (format_sed|format_perl))
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
         if(flags & format_is_if)
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
         if(flags & (format_sed|format_perl))
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
               oi_assign(&out, _reg_format_aux(out, m, fmt, flags | format_is_if, traits_inst));
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*(fmt-1))) == traits_type::syntax_colon)
                  re_skip_format(fmt, traits_inst);
            }
            else
            {
               re_skip_format(fmt, traits_inst);
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*(fmt-1))) == traits_type::syntax_colon)
                  oi_assign(&out, _reg_format_aux(out, m, fmt, flags | format_is_if, traits_inst));
            }
            return out;
         }
      }
      default:
default_opt:
         if((flags & format_sed) && (*fmt == '&'))
         {
            oi_assign(&out, re_copy_out(out, Iterator(m[0].first), Iterator(m[0].second)));
            ++fmt;
            continue;
         }
         *out = *fmt;
         ++out;
         ++fmt;
      }
   }

   return out;
}

#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
} // namespace
#endif

template <class S>
class string_out_iterator
{
   S* out;
public:
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
   unsigned flags;
   const traits_type* pt;

public:
   merge_out_predicate(OutputIterator& o, Iterator& pi, const charT* f, unsigned format_flags, const traits_type& p)
      : out(&o), last(&pi), fmt(f), flags(format_flags), pt(&p){}

   ~merge_out_predicate() {}
   bool BOOST_REGEX_CALL operator()(const boost::match_results<Iterator, Allocator>& m)
   {
      const charT* f = fmt;
      if(0 == (flags & format_no_copy))
         oi_assign(out, re_copy_out(*out, Iterator(m[-1].first), Iterator(m[-1].second)));
      oi_assign(out, _reg_format_aux(*out, m, f, flags, *pt));
      *last = m[-2].first;
      return flags & format_first_only ? false : true;
   }
};

} // namespace re_detail

template <class OutputIterator, class Iterator, class Allocator, class charT>
OutputIterator regex_format(OutputIterator out,
                          const match_results<Iterator, Allocator>& m,
                          const charT* fmt,
                          unsigned flags = 0
                         )
{
   regex_traits<charT> t;
   return re_detail::_reg_format_aux(out, m, fmt, flags, t);
}

template <class OutputIterator, class Iterator, class Allocator, class charT>
OutputIterator regex_format(OutputIterator out,
                          const match_results<Iterator, Allocator>& m,
                          const std::basic_string<charT>& fmt,
                          unsigned flags = 0
                         )
{
   regex_traits<charT> t;
   const charT* start = fmt.c_str();
   return re_detail::_reg_format_aux(out, m, start, flags, t);
}  

template <class Iterator, class Allocator, class charT>
std::basic_string<charT> regex_format(const match_results<Iterator, Allocator>& m, const charT* fmt, unsigned flags = 0)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_format(i, m, fmt, flags);
   return result;
}

template <class Iterator, class Allocator, class charT>
std::basic_string<charT> regex_format(const match_results<Iterator, Allocator>& m, const std::basic_string<charT>& fmt, unsigned flags = 0)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_format(i, m, fmt.c_str(), flags);
   return result;
}

template <class OutputIterator, class Iterator, class traits, class Allocator, class charT>
OutputIterator regex_merge(OutputIterator out,
                         Iterator first,
                         Iterator last,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const charT* fmt, 
                         unsigned int flags = match_default)
{
   Iterator l = first;
   re_detail::merge_out_predicate<OutputIterator, Iterator, charT, Allocator, traits> oi(out, l, fmt, flags, e.get_traits());
   regex_grep(oi, first, last, e, flags);
   return (flags & format_no_copy) ? out : re_detail::re_copy_out(out, l, last);
}

template <class OutputIterator, class Iterator, class traits, class Allocator, class charT>
inline OutputIterator regex_merge(OutputIterator out,
                         Iterator first,
                         Iterator last,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const std::basic_string<charT>& fmt,
                         unsigned int flags = match_default)
{
   return regex_merge(out, first, last, e, fmt.c_str(), flags);
}

template <class traits, class Allocator, class charT>
std::basic_string<charT> regex_merge(const std::basic_string<charT>& s,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const charT* fmt,
                         unsigned int flags = match_default)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_merge(i, s.begin(), s.end(), e, fmt, flags);
   return result;
}

template <class traits, class Allocator, class charT>
std::basic_string<charT> regex_merge(const std::basic_string<charT>& s,
                         const reg_expression<charT, traits, Allocator>& e, 
                         const std::basic_string<charT>& fmt,
                         unsigned int flags = match_default)
{
   std::basic_string<charT> result;
   re_detail::string_out_iterator<std::basic_string<charT> > i(result);
   regex_merge(i, s.begin(), s.end(), e, fmt.c_str(), flags);
   return result;
}

#ifdef __BORLANDC__
  #pragma option pop
#endif

} // namespace boost

#endif  // BOOST_REGEX_FORMAT_HPP






