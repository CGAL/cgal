/*
 *
 * Copyright (c) 2003
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         regex_iterator.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Provides regex_iterator implementation.
  */

#ifndef BOOST_REGEX_V4_REGEX_ITERATOR_HPP
#define BOOST_REGEX_V4_REGEX_ITERATOR_HPP

#include <boost/shared_ptr.hpp>

namespace boost{

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

template <class BidirectionalIterator, 
          class charT,
          class traits,
          class Allocator>
class regex_iterator_implementation 
{
   typedef basic_regex<charT, traits, Allocator> regex_type;

   match_results<BidirectionalIterator> what;  // current match
   BidirectionalIterator                base;  // start of sequence
   BidirectionalIterator                end;   // end of sequence
   const regex_type*                    pre;   // the expression
   match_flag_type                      flags; // flags for matching

public:
   regex_iterator_implementation(const regex_type* p, BidirectionalIterator last, match_flag_type f)
      : base(), end(last), pre(p), flags(f){}
   bool init(BidirectionalIterator first)
   {
      base = first;
      return regex_search(first, end, what, *pre, flags);
   }
   bool compare(const regex_iterator_implementation& that)
   {
      if(this == &that) return true;
      return (pre == that.pre) && (end == that.end) && (flags == that.flags) && (what[0].first == that.what[0].first) && (what[0].second == that.what[0].second);
   }
   const match_results<BidirectionalIterator>& get()
   { return what; }
   bool next()
   {
      if(what.prefix().first != what[0].second)
         flags |= match_prev_avail;
      BidirectionalIterator next_start = what[0].second;
      match_flag_type f(flags);
      if(!what.length())
         f |= regex_constants::match_not_initial_null;
      bool result = regex_search(next_start, end, what, *pre, f);
      if(result)
         what.set_base(base);
      return result;
   }
};

template <class BidirectionalIterator, 
          class charT = BOOST_DEDUCED_TYPENAME re_detail::regex_iterator_traits<BidirectionalIterator>::value_type,
          class traits = regex_traits<charT>,
          class Allocator = BOOST_DEFAULT_ALLOCATOR(charT) >
class regex_iterator 
{
private:
   typedef regex_iterator_implementation<BidirectionalIterator, charT, traits, Allocator> impl;
   typedef shared_ptr<impl> pimpl;
public:
   typedef          basic_regex<charT, traits, Allocator>                   regex_type;
   typedef          match_results<BidirectionalIterator>                    value_type;
   typedef typename re_detail::regex_iterator_traits<BidirectionalIterator>::difference_type 
                                                                            difference_type;
   typedef          const value_type*                                       pointer;
   typedef          const value_type&                                       reference; 
   typedef          std::forward_iterator_tag                               iterator_category;
   
   regex_iterator(){}
   regex_iterator(BidirectionalIterator a, BidirectionalIterator b, 
                  const regex_type& re, 
                  match_flag_type m = match_default)
                  : pdata(new impl(&re, b, m))
   {
      if(!pdata->init(a))
      {
         pdata.reset();
      }
   }
   regex_iterator(const regex_iterator& that)
      : pdata(that.pdata) {}
   regex_iterator& operator=(const regex_iterator& that)
   {
      pdata = that.pdata;
      return *this;
   }
   bool operator==(const regex_iterator& that)const
   { 
      if((pdata.get() == 0) || (that.pdata.get() == 0))
         return pdata.get() == that.pdata.get();
      return pdata->compare(*(that.pdata.get())); 
   }
   bool operator!=(const regex_iterator& that)const
   { return !(*this == that); }
   const value_type& operator*()const
   { return pdata->get(); }
   const value_type* operator->()const
   { return &(pdata->get()); }
   regex_iterator& operator++()
   {
      cow();
      if(0 == pdata->next())
      {
         pdata.reset();
      }
      return *this;
   }
   regex_iterator operator++(int)
   {
      regex_iterator result(*this);
      ++(*this);
      return result;
   }
private:

   pimpl pdata;

   void cow()
   {
      // copy-on-write
      if(pdata.get() && !pdata.unique())
      {
         pdata.reset(new impl(*(pdata.get())));
      }
   }
};

typedef regex_iterator<const char*> cregex_iterator;
typedef regex_iterator<std::string::const_iterator> sregex_iterator;
#ifndef BOOST_NO_WREGEX
typedef regex_iterator<const wchar_t*> wcregex_iterator;
typedef regex_iterator<std::wstring::const_iterator> wsregex_iterator;
#endif

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost

#endif // BOOST_REGEX_V4_REGEX_ITERATOR_HPP

