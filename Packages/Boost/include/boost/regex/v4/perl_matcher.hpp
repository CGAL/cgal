/*
 *
 * Copyright (c) 2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#ifndef BOOST_REGEX_MATCHER_HPP
#define BOOST_REGEX_MATCHER_HPP

#include <boost/regex/v4/iterator_category.hpp>

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost{
namespace re_detail{

//
// error checking API:
//
BOOST_REGEX_DECL void BOOST_REGEX_CALL verify_options(boost::regex::flag_type ef, match_flag_type mf);


//
// Unfortunately Rogue Waves standard library appears to have a bug
// in std::basic_string::compare that results in eroneous answers
// in some cases (tested with Borland C++ 5.1, Rogue Wave lib version
// 0x020101) the test case was:
// {39135,0} < {0xff,0}
// which succeeds when it should not.
//
#ifndef _RWSTD_VER
# define STR_COMP(s,p) s.compare(p)
#else
template <class C, class T, class A>
inline int string_compare(const std::basic_string<C,T,A>& s, const C* p)
{ return s.compare(p); }
inline int string_compare(const std::string& s, const char* p)
{ return std::strcmp(s.c_str(), p); }
# ifndef BOOST_NO_WREGEX
inline int string_compare(const std::wstring& s, const wchar_t* p)
{ return std::wcscmp(s.c_str(), p); }
#endif
# define STR_COMP(s,p) string_compare(s,p)
#endif

template<class charT>
inline const charT* re_skip_past_null(const charT* p)
{
  while (*p != 0) ++p;
  return ++p;
}

template <class iterator, class charT, class traits_type, class Allocator>
iterator BOOST_REGEX_CALL re_is_set_member(iterator next, 
                          iterator last, 
                          const re_set_long* set_, 
                          const reg_expression<charT, traits_type, Allocator>& e)
{   
   const charT* p = reinterpret_cast<const charT*>(set_+1);
   iterator ptr;
   unsigned int i;
   bool icase = e.flags() & regex_constants::icase;

   if(next == last) return next;

   typedef typename traits_type::string_type traits_string_type;
   const traits_type& traits_inst = e.get_traits();
   
   // dwa 9/13/00 suppress incorrect MSVC warning - it claims this is never
   // referenced
   (void)traits_inst;

   // try and match a single character, could be a multi-character
   // collating element...
   for(i = 0; i < set_->csingles; ++i)
   {
      ptr = next;
      if(*p == 0)
      {
         // treat null string as special case:
         if(traits_inst.translate(*ptr, icase) != *p)
         {
            while(*p == 0)++p;
            continue;
         }
         return set_->isnot ? next : (ptr == next) ? ++next : ptr;
      }
      else
      {
         while(*p && (ptr != last))
         {
            if(traits_inst.translate(*ptr, icase) != *p)
               break;
            ++p;
            ++ptr;
         }

         if(*p == 0) // if null we've matched
            return set_->isnot ? next : (ptr == next) ? ++next : ptr;

         p = re_skip_past_null(p);     // skip null
      }
   }

   charT col = traits_inst.translate(*next, icase);


   if(set_->cranges || set_->cequivalents)
   {
      traits_string_type s2(1, col);
      traits_string_type s1;
      //
      // try and match a range, NB only a single character can match
      if(set_->cranges)
      {
         if((e.flags() & regex_constants::collate) == 0)
            s1 = s2;
         else
            traits_inst.transform(s1, s2);
         for(i = 0; i < set_->cranges; ++i)
         {
            if(STR_COMP(s1, p) <= 0)
            {
               while(*p)++p;
               ++p;
               if(STR_COMP(s1, p) >= 0)
                  return set_->isnot ? next : ++next;
            }
            else
            {
               // skip first string
               while(*p)++p;
               ++p;
            }
            // skip second string
            while(*p)++p;
            ++p;
         }
      }
      //
      // try and match an equivalence class, NB only a single character can match
      if(set_->cequivalents)
      {
         traits_inst.transform_primary(s1, s2);
         for(i = 0; i < set_->cequivalents; ++i)
         {
            if(STR_COMP(s1, p) == 0)
               return set_->isnot ? next : ++next;
            // skip string
            while(*p)++p;
            ++p;
         }
      }
   }
   if(traits_inst.is_class(col, set_->cclasses) == true)
      return set_->isnot ? next : ++next;
   return set_->isnot ? ++next : next;
}

template <class charT, class traits, class Allocator>
struct access_t : public reg_expression<charT, traits, Allocator>
{
   typedef typename is_byte<charT>::width_type width_type;
   typedef reg_expression<charT, traits, Allocator> base_type;
   typedef charT char_type;
   typedef traits traits_type;
   typedef Allocator alloc_type;

   static int repeat_count(const base_type& b) 
   { return base_type::repeat_count(b); }
   static unsigned int restart_type(const base_type& b) 
   { return base_type::restart_type(b); }
   static const re_syntax_base* first(const base_type& b)
   { return base_type::first(b); }
   static const unsigned char* get_map(const base_type& b)
   { return base_type::get_map(b); }
   static std::size_t leading_length(const base_type& b)
   { return base_type::leading_length(b); }
   static const kmp_info<charT>* get_kmp(const base_type& b)
   { return base_type::get_kmp(b); }
   static bool can_start(char_type c, const unsigned char* _map, unsigned char mask)
   {
      return reg_expression<char_type, traits_type, alloc_type>::can_start(c, _map, mask, width_type());
   }
};


template <class BidiIterator>
class repeater_count
{
   repeater_count** stack;
   repeater_count* next;
   int id;
   unsigned count;              // the number of iterations so far
   BidiIterator start_pos; // where the last repeat started
public:
   repeater_count(repeater_count** s)
   {
      stack = s;
      next = 0;
      id = -1;
      count = 0;
   }
   repeater_count(int i, repeater_count** s, BidiIterator start)
      : start_pos(start)
   {
      id = i;
      stack = s;
      next = *stack;
      *stack = this;
      if(id > next->id)
         count = 0;
      else
      {
         repeater_count* p = next;
         while(p->id != id)
            p = p->next;
         count = p->count;
         start_pos = p->start_pos;
      }
   }
   ~repeater_count()
   {
      *stack = next;
   }
   unsigned get_count() { return count; }
   int get_id() { return id; }
   int operator++() { return ++count; }
   bool check_null_repeat(const BidiIterator& pos, unsigned max)
   {
      // this is called when we are about to start a new repeat,
      // if the last one was NULL move our count to max,
      // otherwise save the current position.
      bool result = (count == 0) ? false : (pos == start_pos);
      if(result)
         count = max;
      else
         start_pos = pos;
      return result;
   }
};

struct saved_state;

enum saved_state_type
{
   saved_type_end = 0,
   saved_type_paren = 1,
   saved_type_recurse = 2,
   saved_type_assertion = 3,
   saved_state_alt = 4,
   saved_state_repeater_count = 5,
   saved_state_extra_block = 6,
   saved_state_greedy_single_repeat = 7,
   saved_state_rep_slow_dot = 8,
   saved_state_rep_fast_dot = 9,
   saved_state_rep_char = 10,
   saved_state_rep_short_set = 11,
   saved_state_rep_long_set = 12,
   saved_state_non_greedy_long_repeat = 13, 
   saved_state_count = 14
};

template <class BidiIterator, class Allocator, class traits, class Allocator2>
class perl_matcher
{
public:
   typedef typename traits::char_type char_type;
   typedef perl_matcher<BidiIterator, Allocator, traits, Allocator2> self_type;
   typedef bool (self_type::*matcher_proc_type)(void);
   typedef access_t<char_type, traits, Allocator2> access;
   typedef typename traits::size_type traits_size_type;
   typedef typename traits::uchar_type traits_uchar_type;
   typedef typename is_byte<char_type>::width_type width_type;
   typedef typename regex_iterator_traits<BidiIterator>::difference_type difference_type;

   perl_matcher(BidiIterator first, BidiIterator end, 
      match_results<BidiIterator, Allocator>& what, 
      const reg_expression<char_type, traits, Allocator2>& e,
      match_flag_type f);

   bool match();
   bool match_imp();
   bool find();
   bool find_imp();
#ifdef BOOST_REGEX_HAS_MS_STACK_GUARD
   typedef bool (perl_matcher::*protected_proc_type)();
   bool protected_call(protected_proc_type);
#endif

   void setf(match_flag_type f)
   { m_match_flags |= f; }
   void unsetf(match_flag_type f)
   { m_match_flags &= ~f; }

private:
   void estimate_max_state_count(std::random_access_iterator_tag*);
   void estimate_max_state_count(void*);
   bool match_prefix();
   bool match_all_states();

   // match procs, stored in s_match_vtable:
   bool match_startmark();
   bool match_endmark();
   bool match_literal();
   bool match_start_line();
   bool match_end_line();
   bool match_wild();
   bool match_match();
   bool match_word_boundary();
   bool match_within_word();
   bool match_word_start();
   bool match_word_end();
   bool match_buffer_start();
   bool match_buffer_end();
   bool match_backref();
   bool match_long_set();
   bool match_set();
   bool match_jump();
   bool match_alt();
   bool match_rep();
   bool match_combining();
   bool match_soft_buffer_end();
   bool match_restart_continue();
   bool match_long_set_repeat();
   bool match_set_repeat();
   bool match_char_repeat();
   bool match_dot_repeat_fast();
   bool match_dot_repeat_slow();
   bool backtrack_till_match(unsigned count);

   // find procs stored in s_find_vtable:
   bool find_restart_any();
   bool find_restart_word();
   bool find_restart_line();
   bool find_restart_buf();
   bool find_restart_lit();

private:
   // final result structure to be filled in:
   match_results<BidiIterator, Allocator>& m_result;
   // temporary result for POSIX matches:
   scoped_ptr<match_results<BidiIterator, Allocator> > m_temp_match;
   // pointer to actual result structure to fill in:
   match_results<BidiIterator, Allocator>* m_presult;
   // start of sequence being searched:
   BidiIterator base;
   // end of sequence being searched:
   BidiIterator last; 
   // current character being examined:
   BidiIterator position;
   // where to restart next search after failed match attempt:
   BidiIterator restart;
   // where the current search started from, acts as base for $` during grep:
   BidiIterator search_base;
   // the expression being examined:
   const reg_expression<char_type, traits, Allocator2>& re;
   // the expression's traits class:
   const traits& traits_inst;
   // the next state in the machine being matched:
   const re_syntax_base* pstate;
   // matching flags in use:
   match_flag_type m_match_flags;
   // how many states we have examined so far:
   difference_type state_count;
   // max number of states to examine before giving up:
   difference_type max_state_count;
   // whether we should ignore case or not:
   bool icase;
   // set to true when (position == last), indicates that we may have a partial match:
   bool m_has_partial_match;
   // set to true whenever we get a match:
   bool m_has_found_match;
   // the current repeat being examined:
   repeater_count<BidiIterator>* next_count;
   // the first repeat being examined (top of linked list):
   repeater_count<BidiIterator> rep_obj;

#ifdef BOOST_REGEX_NON_RECURSIVE
   //
   // additional members for non-recursive version:
   //
   typedef bool (self_type::*unwind_proc_type)(bool);

   void extend_stack();
   bool unwind(bool);
   bool unwind_end(bool);
   bool unwind_paren(bool);
   bool unwind_recursion_stopper(bool);
   bool unwind_assertion(bool);
   bool unwind_alt(bool);
   bool unwind_repeater_counter(bool);
   bool unwind_extra_block(bool);
   bool unwind_greedy_single_repeat(bool);
   bool unwind_slow_dot_repeat(bool);
   bool unwind_fast_dot_repeat(bool);
   bool unwind_char_repeat(bool);
   bool unwind_short_set_repeat(bool);
   bool unwind_long_set_repeat(bool);
   bool unwind_non_greedy_repeat(bool);
   void destroy_single_repeat();
   void push_matched_paren(int index, const sub_match<BidiIterator>& sub);
   void push_recursion_stopper();
   void push_assertion(const re_syntax_base* ps, bool positive);
   void push_alt(const re_syntax_base* ps);
   void push_repeater_count(int i, repeater_count<BidiIterator>** s);
   void push_single_repeat(unsigned c, const re_repeat* r, BidiIterator last_position, int id);
   void push_non_greedy_repeat(const re_syntax_base* ps);


   // pointer to base of stack:
   saved_state* m_stack_base;
   // pointer to current stack position:
   saved_state* m_backup_state;
   // determines what value to return when unwinding from recursion,
   // allows for mixed recursive/non-recursive algorithm:
   bool m_recursive_result;
   // how many memory blocks have we used up?:
   unsigned used_block_count;
#endif

   // these operations aren't allowed, so are declared private:
   perl_matcher& operator=(const perl_matcher&);
   perl_matcher(const perl_matcher&);
};

} // namespace re_detail

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost

//
// include the implementation of perl_matcher:
//
#ifdef BOOST_REGEX_RECURSIVE
#include <boost/regex/v4/perl_matcher_recursive.hpp>
#else
#include <boost/regex/v4/perl_matcher_non_recursive.hpp>
#endif
// this one has to be last:
#include <boost/regex/v4/perl_matcher_common.hpp>

#endif

