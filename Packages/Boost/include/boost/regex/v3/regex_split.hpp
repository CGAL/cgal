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
  *   FILE         regex_split.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Implements regex_split and associated functions.
  *                Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */

#ifndef BOOST_REGEX_SPLIT_HPP
#define BOOST_REGEX_SPLIT_HPP

namespace boost{

#ifdef __BORLANDC__
   #pragma option push -a8 -b -Vx -Ve -pc
#endif

namespace re_detail{

template <class charT>
const reg_expression<charT>& get_default_expression(charT)
{
   static const charT expression_text[4] = { '\\', 's', '+', '\00', };
   static const reg_expression<charT> e(expression_text);
   return e;
}

template <class OutputIterator, class charT, class Traits1, class Alloc1, class Alloc2>
class split_pred
{
   typedef std::basic_string<charT, Traits1, Alloc1> string_type;
   typedef typename string_type::const_iterator iterator_type;
   iterator_type* p_last;
   OutputIterator* p_out;
   std::size_t* p_max;
   std::size_t initial_max;
public:
   split_pred(iterator_type* a, OutputIterator* b, std::size_t* c)
      : p_last(a), p_out(b), p_max(c), initial_max(*c) {}

   bool operator()(const match_results<iterator_type, Alloc2>& what);
};

template <class OutputIterator, class charT, class Traits1, class Alloc1, class Alloc2>
bool split_pred<OutputIterator, charT, Traits1, Alloc1, Alloc2>::operator()
   (const match_results<iterator_type, Alloc2>& what)
{
   *p_last = what[0].second;
   if(what.size() > 1)
   {
      // output sub-expressions only:
      for(unsigned i = 1; i < what.size(); ++i)
      {
         *(*p_out) = static_cast<string_type>(what[i]);
         ++(*p_out);
         if(0 == --*p_max) return false;
      }
      return *p_max != 0;
   }
   else
   {
      // output $` only if it's not-null or not at the start of the input:
      const sub_match<iterator_type>& sub = what[-1];
      if((sub.first != sub.second) || (*p_max != initial_max))
      {
         *(*p_out) = static_cast<string_type>(sub);
         ++(*p_out);
         return --*p_max;
      }
   }
   //
   // initial null, do nothing:
   return true;
}

} // namespace re_detail

template <class OutputIterator, class charT, class Traits1, class Alloc1, class Traits2, class Alloc2>
std::size_t regex_split(OutputIterator out,
                   std::basic_string<charT, Traits1, Alloc1>& s, 
                   const reg_expression<charT, Traits2, Alloc2>& e,
                   unsigned flags,
                   std::size_t max_split)
{
   typedef typename std::basic_string<charT, Traits1, Alloc1>::const_iterator ci_t;
   ci_t last = s.begin();
   std::size_t init_size = max_split;
   re_detail::split_pred<OutputIterator, charT, Traits1, Alloc1, Alloc2> pred(&last, &out, &max_split);
   ci_t i, j;
   i = s.begin();
   j = s.end();
   regex_grep(pred, i, j, e, flags);
   //
   // if there is still input left, do a final push as long as max_split
   // is not exhausted, and we're not splitting sub-expressions rather 
   // than whitespace:
   if(max_split && (last != s.end()) && (e.mark_count() == 1))
   {
      *out = std::basic_string<charT, Traits1, Alloc1>((ci_t)last, (ci_t)s.end());
      ++out;
      last = s.end();
      --max_split;
   }
   //
   // delete from the string everything that has been processed so far:
   s.erase(0, last - s.begin());
   //
   // return the number of new records pushed:
   return init_size - max_split;
}

template <class OutputIterator, class charT, class Traits1, class Alloc1, class Traits2, class Alloc2>
inline std::size_t regex_split(OutputIterator out,
                   std::basic_string<charT, Traits1, Alloc1>& s, 
                   const reg_expression<charT, Traits2, Alloc2>& e,
                   unsigned flags = match_default)
{
   return regex_split(out, s, e, flags, UINT_MAX);
}

template <class OutputIterator, class charT, class Traits1, class Alloc1>
inline std::size_t regex_split(OutputIterator out,
                   std::basic_string<charT, Traits1, Alloc1>& s)
{
   return regex_split(out, s, re_detail::get_default_expression(charT(0)), match_default, UINT_MAX);
}

#ifdef __BORLANDC__
  #pragma option pop
#endif

} // namespace boost

#endif



