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
  *   FILE         regex_match.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Regular expression matching algorithms.
  *                Note this is an internal header file included
  *                by regex.hpp, do not include on its own.
  */


#ifndef BOOST_REGEX_MATCH_HPP
#define BOOST_REGEX_MATCH_HPP

#ifndef BOOST_REGEX_MAX_STATE_COUNT
#  define BOOST_REGEX_MAX_STATE_COUNT 100000000
#endif

#include <boost/limits.hpp>


namespace boost{
   namespace re_detail{

#ifdef __BORLANDC__
   #pragma option push -a8 -b -Vx -Ve -pc  -w-8026 -w-8027
#endif

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
# endif
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
   bool icase = e.flags() & regbase::icase;

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
         if(e.flags() & regbase::nocollate)
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

template <class iterator, class Allocator>
class _priv_match_data
{
public:
   typedef typename boost::detail::rebind_allocator<int, Allocator>::type i_alloc;
   typedef typename boost::detail::rebind_allocator<iterator, Allocator>::type it_alloc;
   typedef typename regex_iterator_traits<iterator>::difference_type difference_type;

   match_results_base<iterator, Allocator> temp_match;
   // failure stacks:
   jstack<match_results_base<iterator, Allocator>, Allocator> matches;
   jstack<iterator, Allocator> prev_pos;
   jstack<const re_syntax_base*, Allocator> prev_record;
   jstack<int, Allocator> prev_acc;
   int* accumulators;
   unsigned int caccumulators;
   difference_type state_count;
   difference_type max_state_count;
   iterator* loop_starts;

   _priv_match_data(const match_results_base<iterator, Allocator>&, iterator, iterator, std::size_t);
   
   ~_priv_match_data()
   {
      m_free();
   }
   void m_free();
   void set_accumulator_size(unsigned int size);
   int* get_accumulators()
   {
      return accumulators;
   }
   iterator* get_loop_starts()
   {
      return loop_starts;
   }
   void estimate_max_state_count(iterator a, iterator b, std::size_t states, std::random_access_iterator_tag*)
   {
      difference_type dist = boost::re_detail::distance(a,b);
      states *= states;
      difference_type lim = (std::numeric_limits<difference_type>::max)() - 1000 - states;
      if(dist > (difference_type)(lim / states))
         max_state_count = lim;
      else
         max_state_count = 1000 + states * dist;
   }
   void estimate_max_state_count(iterator a, iterator b, std::size_t states, void*)
   {
      // we don't know how long the sequence is:
      max_state_count = BOOST_REGEX_MAX_STATE_COUNT;
   }
};

template <class iterator, class Allocator>
_priv_match_data<iterator, Allocator>::_priv_match_data(const match_results_base<iterator, Allocator>& m, iterator a, iterator b, std::size_t states)
  : temp_match(m), matches(64, m.allocator()), prev_pos(64, m.allocator()), prev_record(64, m.allocator())
{
  typedef typename regex_iterator_traits<iterator>::iterator_category category;
  
  accumulators = 0;
  caccumulators = 0;
  loop_starts = 0;
  state_count = 0;
  estimate_max_state_count(a, b, states, static_cast<category*>(0));
}

template <class iterator, class Allocator>
void _priv_match_data<iterator, Allocator>::set_accumulator_size(unsigned int size)
{
   if(size > caccumulators)
   {
      m_free();
      caccumulators = size;
      accumulators = i_alloc(temp_match.allocator()).allocate(caccumulators);
      BOOST_REGEX_NOEH_ASSERT(accumulators)
      loop_starts = it_alloc(temp_match.allocator()).allocate(caccumulators);
      BOOST_REGEX_NOEH_ASSERT(loop_starts)
      for(unsigned i = 0; i < caccumulators; ++i)
         new (loop_starts + i) iterator();
   }
}

template <class iterator, class Allocator>
void _priv_match_data<iterator, Allocator>::m_free()
{
   if(caccumulators)
   {
      i_alloc temp1(temp_match.allocator());
      temp1.deallocate(accumulators, caccumulators);
      for(unsigned i = 0; i < caccumulators; ++i)
         ::boost::re_detail::pointer_destroy(loop_starts + i);
      it_alloc temp2(temp_match.allocator());
      temp2.deallocate(loop_starts, caccumulators);
   }
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


#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
//
// Ugly ugly hack,
// template don't merge if they contain switch statements so declare these
// templates in unnamed namespace (ie with internal linkage), each translation
// unit then gets its own local copy, it works seemlessly but bloats the app.
namespace{
#endif

template <class iterator, class Allocator, class charT, class traits, class Allocator2>
bool query_match_aux(iterator first, 
                     iterator last, 
                     match_results<iterator, Allocator>& m, 
                     const reg_expression<charT, traits, Allocator2>& e, 
                     unsigned flags,
                     _priv_match_data<iterator, Allocator>& pd,
                     iterator* restart)
{
   typedef access_t<charT, traits, Allocator2> access;

   if(e.flags() & regbase::failbit)
      return false;

   typedef typename traits::size_type traits_size_type;
   typedef typename traits::uchar_type traits_uchar_type;
   typedef typename is_byte<charT>::width_type width_type;
   typedef typename re_detail::regex_iterator_traits<iterator>::difference_type difference_type;

   // declare some local aliases to reduce pointer loads
   // good optimising compilers should make this unnecessary!!
   jstack<match_results_base<iterator, Allocator>, Allocator>& matches = pd.matches;
   jstack<iterator, Allocator>& prev_pos = pd.prev_pos;
   jstack<const re_syntax_base*, Allocator>& prev_record = pd.prev_record;
   jstack<int, Allocator>& prev_acc = pd.prev_acc;
   match_results_base<iterator, Allocator>& temp_match = pd.temp_match;
   temp_match.set_first(first);
   difference_type& state_count = pd.state_count;

   const re_syntax_base* ptr = access::first(e);
   bool match_found = false;
   bool have_partial_match = false;
   bool unwind_stack = false;
   bool need_push_match = (e.mark_count() > 1);
   int cur_acc = -1;    // no active accumulator
   pd.set_accumulator_size(access::repeat_count(e));
   int* accumulators = pd.get_accumulators();
   iterator* start_loop = pd.get_loop_starts();
   int k; // for loops
   bool icase = e.flags() & regbase::icase;
   *restart = first;
   iterator base = first;
   const traits& traits_inst = e.get_traits();
   // dwa 9/13/00 suppress incorrect MSVC warning - it claims this is never
   // referenced
   (void)traits_inst;

   // prepare m for failure:
   /*
   if((flags & match_init) == 0)
   {
      m.init_fail(first, last);
   } */

   retry:

   while(first != last)
   {
      jm_assert(ptr);
      ++state_count;
      switch(ptr->type)
      {
      case syntax_element_match:
         match_jump:
         {
            // match found, save then fallback in case we missed a
            // longer one.
            if((flags & match_not_null) && (first == temp_match[0].first))
               goto failure;
            if((flags & match_all) && (first != last))
               goto failure;
            temp_match.set_second(first);
            m.maybe_assign(temp_match);
            match_found = true;
            if(((flags & match_any) && ((first == last) || !(flags & match_all))) || ((first == last) && (need_push_match == false)))
            {
               // either we don't care what we match or we've matched
               // the whole string and can't match anything longer.
               while(matches.empty() == false)
                  matches.pop();
               while(prev_pos.empty() == false)
                  prev_pos.pop();
               while(prev_record.empty() == false)
                  prev_record.pop();
               while(prev_acc.empty() == false)
                  prev_acc.pop();
               return true;
            }
         }
         goto failure;
      case syntax_element_startmark:
         start_mark_jump:
         if(static_cast<const re_brace*>(ptr)->index > 0)
         {
            temp_match.set_first(first, static_cast<const re_brace*>(ptr)->index);
         }
         else if(
            (static_cast<const re_brace*>(ptr)->index == -1)
            || (static_cast<const re_brace*>(ptr)->index == -2)
         )
         {
           matches.push(temp_match);
            for(k = 0; k <= cur_acc; ++k)
               prev_pos.push(start_loop[k]);
            prev_pos.push(first);
            prev_record.push(ptr);
            for(k = 0; k <= cur_acc; ++k)
               prev_acc.push(accumulators[k]);
            prev_acc.push(cur_acc);
            prev_acc.push(match_found);
            match_found = false;
            // skip next jump and fall through:
            ptr = ptr->next.p;
         }
         ptr = ptr->next.p;
         break;
      case syntax_element_endmark:
         end_mark_jump:
         if(static_cast<const re_brace*>(ptr)->index > 0)
         {
            temp_match.set_second(first, static_cast<const re_brace*>(ptr)->index);
         }
         else if(
            (static_cast<const re_brace*>(ptr)->index == -1)
            || (static_cast<const re_brace*>(ptr)->index == -2)
         )
         {
            match_found = true;
            unwind_stack = true;
            goto failure;
         }
         ptr = ptr->next.p;
         break;
      case syntax_element_literal:
      {
         unsigned int len = static_cast<const re_literal*>(ptr)->length;
         const charT* what = reinterpret_cast<const charT*>(static_cast<const re_literal*>(ptr) + 1);
         //
         // compare string with what we stored in
         // our records:
         for(unsigned int i = 0; i < len; ++i, ++first)
         {
            if((first == last) || (traits_inst.translate(*first, icase) != what[i]))
               goto failure;
         }
         ptr = ptr->next.p;
         break;
      }
      case syntax_element_start_line:
         outer_line_check:
         if(first == temp_match[0].first)
         {
            // we're at the start of the buffer
            if(flags & match_prev_avail)
            {
               inner_line_check:
               // check the previous value even though its before
               // the start of our "buffer".
               iterator t(first);
               --t;
               if(traits::is_separator(*t) && !((*t == '\r') && (*first == '\n')) )
               {
                  ptr = ptr->next.p;
                  continue;
               }
               goto failure;
            }
            if((flags & match_not_bol) == 0)
            {
               ptr = ptr->next.p;
               continue;
            }
            goto failure;
         }
         // we're in the middle of the string
         goto inner_line_check;
      case syntax_element_end_line:
         // we're not yet at the end so *first is always valid:
         if(traits::is_separator(*first))
         {
            if((first != base) || (flags & match_prev_avail))
            {
               // check that we're not in the middle of \r\n sequence
               iterator t(first);
               --t;
               if((*t == '\r') && (*first == '\n'))
               {
                  goto failure;
               }
            }
            ptr = ptr->next.p;
            continue;
         }
         goto failure;
      case syntax_element_wild:
         // anything except possibly NULL or \n:
         if(traits::is_separator(*first))
         {
            if(flags & match_not_dot_newline)
               goto failure;
            ptr = ptr->next.p;
            ++first;
            continue;
         }
         if(*first == charT(0))
         {
            if(flags & match_not_dot_null)
               goto failure;
            ptr = ptr->next.p;
            ++first;
            continue;
         }
         ptr = ptr->next.p;
         ++first;
         break;
      case syntax_element_word_boundary:
      {
         // prev and this character must be opposites:
#if defined(BOOST_REGEX_USE_C_LOCALE) && defined(__GNUC__) && (__GNUC__ == 2) && (__GNUC_MINOR__ < 95)
         bool b = traits::is_class(*first, traits::char_class_word);
#else
         bool b = traits_inst.is_class(*first, traits::char_class_word);
#endif
         if((first == temp_match[0].first)  && ((flags & match_prev_avail) == 0))
         {
            if(flags & match_not_bow)
               b ^= true;
            else
               b ^= false;
         }
         else
         {
            --first;
            b ^= traits_inst.is_class(*first, traits::char_class_word);
            ++first;
         }
         if(b)
         {
            ptr = ptr->next.p;
            continue;
         }
         goto failure;
      }
      case syntax_element_within_word:
         // both prev and this character must be traits::char_class_word:
         if(traits_inst.is_class(*first, traits::char_class_word))
         {
            bool b;
            if((first == temp_match[0].first) && ((flags & match_prev_avail) == 0))
               b = false;
            else
            {
               --first;
               b = traits_inst.is_class(*first, traits::char_class_word);
               ++first;
            }
            if(b)
            {
               ptr = ptr->next.p;
               continue;
            }
         }
         goto failure;
      case syntax_element_word_start:
         if((first == temp_match[0].first) && ((flags & match_prev_avail) == 0))
         {
            // start of buffer:
            if(flags & match_not_bow)
               goto failure;
            if(traits_inst.is_class(*first, traits::char_class_word))
            {
               ptr = ptr->next.p;
               continue;
            }
            goto failure;
         }
         // otherwise inside buffer:
         if(traits_inst.is_class(*first, traits::char_class_word))
         {
            iterator t(first);
            --t;
            if(traits_inst.is_class(*t, traits::char_class_word) == false)
            {
               ptr = ptr->next.p;
               continue;
            }
         }
         goto failure;      // if we fall through to here then we've failed
      case syntax_element_word_end:
         if((first == temp_match[0].first) && ((flags & match_prev_avail) == 0))
            goto failure;  // start of buffer can't be end of word

         // otherwise inside buffer:
         if(traits_inst.is_class(*first, traits::char_class_word) == false)
         {
            iterator t(first);
            --t;
            if(traits_inst.is_class(*t, traits::char_class_word))
            {
               ptr = ptr->next.p;
               continue;
            }
         }
         goto failure;      // if we fall through to here then we've failed
      case syntax_element_buffer_start:
         if((first != temp_match[0].first) || (flags & match_not_bob))
            goto failure;
         // OK match:
         ptr = ptr->next.p;
         break;
      case syntax_element_buffer_end:
         if((first != last) || (flags & match_not_eob))
            goto failure;
         // OK match:
         ptr = ptr->next.p;
         break;
      case syntax_element_backref:
      {
         // compare with what we previously matched:
         iterator i = temp_match[static_cast<const re_brace*>(ptr)->index].first;
         iterator j = temp_match[static_cast<const re_brace*>(ptr)->index].second;
         while(i != j)
         {
            if((first == last) || (traits_inst.translate(*first, icase) != traits_inst.translate(*i, icase)))
               goto failure;
            ++i;
            ++first;
         }
         ptr = ptr->next.p;
         break;
      }
      case syntax_element_long_set:
      {
         // let the traits class do the work:
         iterator t = re_is_set_member(first, last, static_cast<const re_set_long*>(ptr), e);
         if(t != first)
         {
            ptr = ptr->next.p;
            first = t;
            continue;
         }
         goto failure;
      }
      case syntax_element_set:
         // lookup character in table:
         if(static_cast<const re_set*>(ptr)->_map[(traits_uchar_type)traits_inst.translate(*first, icase)])
         {
            ptr = ptr->next.p;
            ++first;
            continue;
         }
         goto failure;
      case syntax_element_jump:
         ptr = static_cast<const re_jump*>(ptr)->alt.p;
         continue;
      case syntax_element_alt:
      {
         // alt_jump:
         if(access::can_start(*first, static_cast<const re_jump*>(ptr)->_map, (unsigned char)mask_take))
         {
            // we can take the first alternative,
            // see if we need to push next alternative:
            if(access::can_start(*first, static_cast<const re_jump*>(ptr)->_map, mask_skip))
            {
               if(need_push_match)
                  matches.push(temp_match);
               for(k = 0; k <= cur_acc; ++k)
                  prev_pos.push(start_loop[k]);
               prev_pos.push(first);
               prev_record.push(ptr);
               for(k = 0; k <= cur_acc; ++k)
                  prev_acc.push(accumulators[k]);
               prev_acc.push(cur_acc);
            }
            ptr = ptr->next.p;
            continue;
         }
         if(access::can_start(*first, static_cast<const re_jump*>(ptr)->_map, mask_skip))
         {
            ptr = static_cast<const re_jump*>(ptr)->alt.p;
            continue;
         }
         goto failure;  // neither option is possible
      }
      case syntax_element_rep:
      {
         // repeater_jump:
         // if we're moving to a higher id (nested repeats etc)
         // zero out our accumualtors:
         if(cur_acc < static_cast<const re_repeat*>(ptr)->id)
         {
            cur_acc = static_cast<const re_repeat*>(ptr)->id;
            accumulators[cur_acc] = 0;
            start_loop[cur_acc] = first;
         }

         cur_acc = static_cast<const re_repeat*>(ptr)->id;

         if(static_cast<const re_repeat*>(ptr)->leading)
            *restart = first;

         //charT c = traits_inst.translate(*first);

         // first of all test for special case where this is last element,
         // if that is the case then repeat as many times as possible,
         // as long as the repeat is greedy:

         if((static_cast<const re_repeat*>(ptr)->alt.p->type == syntax_element_match)
            && (static_cast<const re_repeat*>(ptr)->greedy == true))
         {
            // see if we can take the repeat:
            if(((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->max)
                  && access::can_start(*first, static_cast<const re_repeat*>(ptr)->_map, mask_take))
            {
               // push terminating match as fallback:
               if((unsigned int)accumulators[cur_acc] >= static_cast<const re_repeat*>(ptr)->min)
               {
                  if((prev_record.empty() == false) && (prev_record.peek() == static_cast<const re_repeat*>(ptr)->alt.p))
                  {
                     // we already have the required fallback
                     // don't add any more, just update this one:
                     if(need_push_match)
                        matches.peek() = temp_match;
                     prev_pos.peek() = first;
                  }
                  else
                  {
                     if(need_push_match)
                        matches.push(temp_match);
                     prev_pos.push(first);
                     prev_record.push(static_cast<const re_repeat*>(ptr)->alt.p);
                  }
               }
               // move to next item in list:
               if((first != start_loop[cur_acc]) || !accumulators[cur_acc])
               {
                  ++accumulators[cur_acc];
                  ptr = ptr->next.p;
                  start_loop[cur_acc] = first;
                  continue;
               }
               else if((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->min)
               {
                  // the repeat was null, and we haven't gone round min times yet,
                  // since all subsequent repeats will be null as well, just update
                  // our repeat count and skip out.
                  accumulators[cur_acc] = static_cast<const re_repeat*>(ptr)->min;
                  ptr = static_cast<const re_repeat*>(ptr)->alt.p;
                  continue;
               }
               goto failure;
            }
            // see if we can skip the repeat:
            if(((unsigned int)accumulators[cur_acc] >= static_cast<const re_repeat*>(ptr)->min)
               && access::can_start(*first, static_cast<const re_repeat*>(ptr)->_map, mask_skip))
            {
               ptr = static_cast<const re_repeat*>(ptr)->alt.p;
               continue;
            }
            // otherwise fail:
            goto failure;
         }

         // OK if we get to here then the repeat is either non-terminal or non-greedy,
         // see if we can skip the repeat:
         if(((unsigned int)accumulators[cur_acc] >= static_cast<const re_repeat*>(ptr)->min)
            && access::can_start(*first, static_cast<const re_repeat*>(ptr)->_map, mask_skip))
         {
            // see if we can push failure info:
            if(((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->max)
               && access::can_start(*first, static_cast<const re_repeat*>(ptr)->_map, mask_take))
            {
               // check to see if the last loop matched a NULL string
               // if so then we really don't want to loop again:
               if(((unsigned int)accumulators[cur_acc] == static_cast<const re_repeat*>(ptr)->min)
                  || (first != start_loop[cur_acc]))
               {
                  if(need_push_match)
                     matches.push(temp_match);
                  prev_pos.push(first);
                  prev_record.push(ptr);
                  for(k = 0; k <= cur_acc; ++k)
                     prev_acc.push(accumulators[k]);
                  // for non-greedy repeats save whether we have a match already:
                  if(static_cast<const re_repeat*>(ptr)->greedy == false)
                  {
                     prev_acc.push(match_found);
                     match_found = false;
                  }
               }
            }
            ptr = static_cast<const re_repeat*>(ptr)->alt.p;
            continue;
         }

         // otherwise see if we can take the repeat:
         if(((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->max)
               && access::can_start(*first, static_cast<const re_repeat*>(ptr)->_map, mask_take) &&
               ((first != start_loop[cur_acc]) || !accumulators[cur_acc]))
         {
            // move to next item in list:
            ++accumulators[cur_acc];
            ptr = ptr->next.p;
            start_loop[cur_acc] = first;
            continue;
         }
         else if((first == start_loop[cur_acc]) && accumulators[cur_acc] && ((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->min))
         {
            // the repeat was null, and we haven't gone round min times yet,
            // since all subsequent repeats will be null as well, just update
            // our repeat count and skip out.
            accumulators[cur_acc] = static_cast<const re_repeat*>(ptr)->min;
            ptr = static_cast<const re_repeat*>(ptr)->alt.p;
            continue;
         }

         // if we get here then neither option is allowed so fail:
         goto failure;

      }
      case syntax_element_combining:
         if(traits_inst.is_combining(traits_inst.translate(*first, icase)))
            goto failure;
         ++first;
         while((first != last) && traits_inst.is_combining(traits_inst.translate(*first, icase)))++first;
         ptr = ptr->next.p;
         continue;
      case syntax_element_soft_buffer_end:
         {
            if(flags & match_not_eob)
               goto failure;
            iterator p(first);
            while((p != last) && traits_inst.is_separator(traits_inst.translate(*p, icase)))++p;
            if(p != last)
               goto failure;
            ptr = ptr->next.p;
            continue;
         }
      case syntax_element_restart_continue:
         if(first != temp_match[-1].first)
            goto failure;
         ptr = ptr->next.p;
         continue;
      default:
         jm_assert(0); // should never get to here!!
         return false;
      }
   }

   //
   // if we get to here then we've run out of characters to match against,
   // we could however still have non-character regex items left
   if((ptr->can_be_null == 0) && ((flags & match_partial) == 0))
      goto failure;
   while(true)
   {
      jm_assert(ptr);
      ++state_count;
      switch(ptr->type)
      {
      case syntax_element_match:
         goto match_jump;
      case syntax_element_startmark:
         goto start_mark_jump;
      case syntax_element_endmark:
         goto end_mark_jump;
      case syntax_element_start_line:
         goto outer_line_check;
      case syntax_element_end_line:
         // we're at the end so *first is never valid:
         if((flags & match_not_eol) == 0)
         {
            ptr = ptr->next.p;
            continue;
         }
         goto failure;
      case syntax_element_word_boundary:
      case syntax_element_word_end:
         if(((flags & match_not_eow) == 0) && (first != temp_match[0].first))
         {
            iterator t(first);
            --t;
            if(traits_inst.is_class(*t, traits::char_class_word))
            {
               ptr = ptr->next.p;
               continue;
            }
         }
         goto failure;
      case syntax_element_buffer_end:
      case syntax_element_soft_buffer_end:
         if(flags & match_not_eob)
            goto failure;
         // OK match:
         ptr = ptr->next.p;
         break;
      case syntax_element_jump:
         ptr = static_cast<const re_jump*>(ptr)->alt.p;
         continue;
      case syntax_element_alt:
         if(ptr->can_be_null & mask_take)
         {
            // we can test the first alternative,
            // see if we need to push next alternative:
            if(ptr->can_be_null & mask_skip)
            {
               if(need_push_match)
                  matches.push(temp_match);
               for(k = 0; k <= cur_acc; ++k)
                  prev_pos.push(start_loop[k]);
               prev_pos.push(first);
               prev_record.push(ptr);
               for(k = 0; k <= cur_acc; ++k)
                  prev_acc.push(accumulators[k]);
               prev_acc.push(cur_acc);
            }
            ptr = ptr->next.p;
            continue;
         }
         if(ptr->can_be_null & mask_skip)
         {
            ptr = static_cast<const re_jump*>(ptr)->alt.p;
            continue;
         }
         goto failure;  // neither option is possible
      case syntax_element_rep:
         // if we're moving to a higher id (nested repeats etc)
         // zero out our accumualtors:
         if(cur_acc < static_cast<const re_repeat*>(ptr)->id)
         {
            cur_acc = static_cast<const re_repeat*>(ptr)->id;
            accumulators[cur_acc] = 0;
            start_loop[cur_acc] = first;
         }

         cur_acc = static_cast<const re_repeat*>(ptr)->id;

         // see if we can skip the repeat:
         if(((unsigned int)accumulators[cur_acc] >= static_cast<const re_repeat*>(ptr)->min)
            && ((ptr->can_be_null & mask_skip) || (flags & match_partial)))
         {
            // don't push failure info, there's no point:
            ptr = static_cast<const re_repeat*>(ptr)->alt.p;
            continue;
         }

         // otherwise see if we can take the repeat:
         if(((unsigned int)accumulators[cur_acc] < static_cast<const re_repeat*>(ptr)->max)
               && (((ptr->can_be_null & (mask_take | mask_skip)) == (mask_take | mask_skip))) || (flags & match_partial))
         {
            // move to next item in list:
            ++accumulators[cur_acc];
            ptr = ptr->next.p;
            start_loop[cur_acc] = first;
            continue;
         }

         // if we get here then neither option is allowed so fail:
         goto failure;
      case syntax_element_restart_continue:
         if(first != temp_match[-1].first)
            goto failure;
         ptr = ptr->next.p;
         continue;
      case syntax_element_backref:
         if(temp_match[static_cast<const re_brace*>(ptr)->index].first
               != temp_match[static_cast<const re_brace*>(ptr)->index].second)
               goto failure;
         ptr = ptr->next.p;
         continue;
      default:
         goto failure;
      }
   }

   failure:

   //
   // check to see if we've been searching too many states:
   //
   if(state_count >= pd.max_state_count)
   {
#ifndef BOOST_NO_EXCEPTIONS
      throw std::runtime_error("Max regex search depth exceeded.");
#else
      while(matches.empty() == false)
         matches.pop();
      while(prev_pos.empty() == false)
         prev_pos.pop();
      while(prev_record.empty() == false)
         prev_record.pop();
      while(prev_acc.empty() == false)
         prev_acc.pop();
      return false;
#endif
   }

   //
   // check for possible partial match:
   //
   if((flags & match_partial)
      && !match_found          // no full match already
      && (base != first)       // some charcters have been consumed
      && (first == last))      // end of input has been reached
   {
      have_partial_match = true;
      temp_match.set_second(first, 0, false);
      m.maybe_assign(temp_match);
   }

   if(prev_record.empty() == false)
   {
      ptr = prev_record.peek();
      switch(ptr->type)
      {
      case syntax_element_alt:
         // get next alternative:
         ptr = static_cast<const re_jump*>(ptr)->alt.p;
         if(need_push_match)
            matches.pop(temp_match);
         prev_acc.pop(cur_acc);
         for(k = cur_acc; k >= 0; --k)
            prev_acc.pop(accumulators[k]);
         prev_pos.pop(first);
         for(k = cur_acc; k >= 0; --k)
            prev_pos.pop(start_loop[k]);
         prev_record.pop();
         if(unwind_stack) goto failure; // unwinding forward assert
         goto retry;
      case syntax_element_rep:
      {
         // we're doing least number of repeats first,
         // increment count and repeat again:
         bool saved_matched = match_found;
         if(need_push_match)
            matches.pop(temp_match);
         prev_pos.pop(first);
         cur_acc = static_cast<const re_repeat*>(ptr)->id;
         if(static_cast<const re_repeat*>(ptr)->greedy == false)
         {
            saved_matched = prev_acc.peek();
            prev_acc.pop();
         }
         for(k = cur_acc; k >= 0; --k)
            prev_acc.pop(accumulators[k]);
         prev_record.pop();
         if(unwind_stack) goto failure; // unwinding forward assert
         if((unsigned int)++accumulators[cur_acc] > static_cast<const re_repeat*>(ptr)->max)
            goto failure;  // repetions exhausted.
         //
         // if the repeat is non-greedy, and we found a match then fail again:
         if((static_cast<const re_repeat*>(ptr)->greedy == false) && (match_found == true))
         {
            goto failure;
         }
         else if (match_found == false)
            match_found = saved_matched;
         ptr = ptr->next.p;
         start_loop[cur_acc] = first;
         goto retry;
      }
      case syntax_element_startmark:
      {
         bool saved_matched = match_found;
         matches.pop(temp_match);
         match_found = prev_acc.peek();
         prev_acc.pop();
         prev_acc.pop(cur_acc);
         for(k = cur_acc; k >= 0; --k)
            prev_acc.pop(accumulators[k]);
         prev_pos.pop(first);
         for(k = cur_acc; k >= 0; --k)
            prev_pos.pop(start_loop[k]);
         prev_record.pop();
         unwind_stack = false;
         if(static_cast<const re_brace*>(ptr)->index == -1)
         {
            if (saved_matched == false)
               goto failure;
            ptr = static_cast<const re_jump*>(ptr->next.p)->alt.p->next.p;
            goto retry;
         }
         if(static_cast<const re_brace*>(ptr)->index == -2)
         {
            if (saved_matched == true)
               goto failure;
            ptr = static_cast<const re_jump*>(ptr->next.p)->alt.p->next.p;
            goto retry;
         }
         else goto failure;
      }
      case syntax_element_match:
         if(need_push_match)
            matches.pop(temp_match);
         prev_pos.pop(first);
         prev_record.pop();
         if(unwind_stack) goto failure; // unwinding forward assert
         goto retry;
     default:
         jm_assert(0);
         // mustn't get here!!
      }
   }

   if(match_found || have_partial_match)
   {
      pd.state_count = 0;
      return true;
   }

   // if we get to here then everything has failed
   // and no match was found:
   return false;
}
#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
} // namespace
#endif


template <class iterator>
void _skip_and_inc(unsigned int& clines, iterator& last_line, iterator& first, const iterator last)
{
   while(first != last)
   {
      if(*first == '\n')
      {
         last_line = ++first;
         ++clines;
      }
      else
         ++first;
   }
}

template <class iterator>
void _skip_and_dec(unsigned int& clines, iterator& last_line, iterator& first, iterator base, std::size_t len)
{
   bool need_line = false;
   for(std::size_t i = 0; i < len; ++i)
   {
      --first;
      if(*first == '\n')
      {
         need_line = true;
         --clines;
      }
   }

   if(need_line)
   {
      last_line = first;

      if(last_line != base)
         --last_line;
      else
         return;

      while((last_line != base) && (*last_line != '\n'))
         --last_line;
      if(*last_line == '\n')
         ++last_line;
   }
}

template <class iterator>
inline void _inc_one(unsigned int& clines, iterator& last_line, iterator& first)
{
   if(*first == '\n')
   {
      last_line = ++first;
      ++clines;
   }
   else
      ++first;
}

template <class iterator, class Allocator>
struct grep_search_predicate
{
   match_results<iterator, Allocator>* pm;
   grep_search_predicate(match_results<iterator, Allocator>* p) : pm(p) {}
   bool operator()(const match_results<iterator, Allocator>& m) 
   {
      *pm = static_cast<const match_results_base<iterator, Allocator>&>(m);
      return false;
   }
};

#if !defined(BOOST_NO_EXPLICIT_FUNCTION_TEMPLATE_ARGUMENTS) && !defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)

template <class iterator, class Allocator>
inline const match_results_base<iterator, Allocator>& grep_out_type(const grep_search_predicate<iterator, Allocator>& o, const Allocator&)
{
   return *(o.pm);
}

#endif

template <class T, class Allocator>
inline const Allocator& grep_out_type(const T&, const Allocator& a)
{
   return a;
}

#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
//
// Ugly ugly hack,
// template don't merge if they contain switch statements so declare these
// templates in unnamed namespace (ie with internal linkage), each translation
// unit then gets its own local copy, it works seemlessly but bloats the app.
namespace{
#endif

//
// reg_grep2:
// find all non-overlapping matches within the sequence first last:
//
template <class Predicate, class I, class charT, class traits, class A, class A2>
unsigned int reg_grep2(Predicate foo, I first, I last, const reg_expression<charT, traits, A>& e, unsigned flags, A2 a)
{
   typedef access_t<charT, traits, A> access;

   if(e.flags() & regbase::failbit)
      return 0;

   typedef typename traits::size_type traits_size_type;
   typedef typename traits::uchar_type traits_uchar_type;
   typedef typename is_byte<charT>::width_type width_type;

   match_results<I, A2> m(grep_out_type(foo, a));
   I restart;
   m.set_size(e.mark_count(), first, last);
   m.set_line(1, first);
   m.set_base(first);

   unsigned int clines = 1;
   unsigned int cmatches = 0;
   I last_line = first;
   I next_base;
   I base = first;
   bool need_init;
   bool leading_match = false;
   const traits& traits_inst = e.get_traits();
   // dwa 9/13/00 suppress incorrect MSVC warning - it claims this is never
   // referenced
   (void)traits_inst;

   flags |= match_init;

   _priv_match_data<I, A2> pd(m, first, last, e.size());

   const unsigned char* _map = access::get_map(e);
   unsigned int type;

   if(first == last)
   {
      // special case, only test if can_be_null,
      // don't dereference any pointers!!
      if(access::first(e)->can_be_null)
      {
         if(query_match_aux(first, last, m, e, flags, pd, &restart))
         {
            foo(m);
            ++cmatches;
         }
      }
      return cmatches;
   }

   // try one time whatever:
   if( access::can_start(*first, _map, (unsigned char)mask_any) )
   {
      if(query_match_aux(first, last, m, e, flags, pd, &restart))
      {
         ++cmatches;
         leading_match = true;
         if(foo(m) == false)
            return cmatches;
         if(m[0].second == last)
            return cmatches;
         // update to end of what matched
         // trying to match again with match_not_null set if this 
         // is a null match...
         need_init = true;
         if(first == m[0].second)
         {
            next_base = m[0].second;
            pd.temp_match.init_fail(next_base, last);
            m.init_fail(next_base, last);
            if(query_match_aux(first, last, m, e, flags | match_not_null, pd, &restart))
            {
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
            }
            else
            {
               need_init = false;
               leading_match = false;
               for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                   {} // dwa 10/20/2000 - warning suppression for MWCW
               if(restart != last)
                  ++restart;
               _skip_and_inc(clines, last_line, first, restart);
            }
         }
         if(need_init)
         {
            _skip_and_inc(clines, last_line, first, m[0].second);
            next_base = m[0].second;
            pd.temp_match.init_fail(next_base, last);
            m.init_fail(next_base, last);
         }
      }
      else
      {
         for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
             {} // dwa 10/20/2000 - warning suppression for MWCW
         if(restart != last)
            ++restart;
         _skip_and_inc(clines, last_line, first, restart);
      }
   }
   else
      _inc_one(clines, last_line, first); 
   flags |= match_prev_avail | match_not_bob;

   
   // depending on what the first record is we may be able to
   // optimise the search:
   type = (flags & match_continuous) ? 
      static_cast<unsigned int>(regbase::restart_continue) 
         : static_cast<unsigned int>(access::restart_type(e));

   if(type == regbase::restart_buf)
      return cmatches;

   switch(type)
   {
   case regbase::restart_lit: 
   case regbase::restart_fixed_lit:
   {
      const kmp_info<charT>* info = access::get_kmp(e);
      int len = info->len;
      const charT* x = info->pstr;
      int j = 0; 
      bool icase = e.flags() & regbase::icase;
      while (first != last) 
      {
         while((j > -1) && (x[j] != traits_inst.translate(*first, icase))) 
            j = info->kmp_next[j];
         _inc_one(clines, last_line, first);
         ++j;
         if(j >= len) 
         {
            if(type == regbase::restart_fixed_lit)
            {
               _skip_and_dec(clines, last_line, first, base, j);
               restart = first;
               std::advance(restart, len);
               m.set_first(first);
               m.set_second(restart);
               m.set_line(clines, last_line);
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
               if(m[0].second == last)
                  return cmatches;
               _skip_and_inc(clines, last_line, first, restart);
               next_base = m[0].second;
               pd.temp_match.init_fail(next_base, last);
               m.init_fail(next_base, last);
               j = 0;
            }
            else
            {
               restart = first;
               _skip_and_dec(clines, last_line, first, base, j);
               if(query_match_aux(first, last, m, e, flags, pd, &restart))
               {

                  m.set_line(clines, last_line);
                  ++cmatches;
                  if(foo(m) == false)
                     return cmatches;
                  if(m[0].second == last)
                     return cmatches;
                  // update to end of what matched
                 _skip_and_inc(clines, last_line, first, m[0].second);
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
                  j = 0;
               }
               else
               {
                  for(int k = 0; (restart != first) && (k < j); ++k, --restart)
                      {} // dwa 10/20/2000 - warning suppression for MWCW
                  if(restart != last)
                     ++restart;
                  _skip_and_inc(clines, last_line, first, restart);
                  j = 0;  //we could do better than this...
               }
            }
         }
      }
      break;
   }
   case regbase::restart_any:
   {
      while(first != last)
      {
         if( access::can_start(*first, _map, (unsigned char)mask_any) )
         {
            if(query_match_aux(first, last, m, e, flags, pd, &restart))
            {

               m.set_line(clines, last_line);
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
               if(m[0].second == last)
                  return cmatches;
               // update to end of what matched
               // trying to match again with match_not_null set if this 
               // is a null match...
               need_init = true;
               if(first == m[0].second)
               {
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
                  if(query_match_aux(first, last, m, e, flags | match_not_null, pd, &restart))
                  {
                     m.set_line(clines, last_line);
                     ++cmatches;
                     if(foo(m) == false)
                        return cmatches;
                  }
                  else
                  {
                     need_init = false;
                     for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                         {} // dwa 10/20/2000 - warning suppression for MWCW
                     if(restart != last)
                        ++restart;
                     _skip_and_inc(clines, last_line, first, restart);
                  }
               }
               if(need_init)
               {
                 _skip_and_inc(clines, last_line, first, m[0].second);
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
               }
               continue;
            }
            else
            {
               for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                   {} // dwa 10/20/2000 - warning suppression for MWCW
               if(restart != last)
                  ++restart;
               _skip_and_inc(clines, last_line, first, restart);
            }
         }
         else
            _inc_one(clines, last_line, first);
      }
   }
   break;
   case regbase::restart_word:
   {
      // do search optimised for word starts:
      while(first != last)
      {
         --first;
         if(*first == '\n')
            --clines;
         // skip the word characters:
         while((first != last) && traits_inst.is_class(*first, traits::char_class_word))
            ++first;
         // now skip the white space:
         while((first != last) && (traits_inst.is_class(*first, traits::char_class_word) == false))
         {
         #ifdef __GNUC__
            //
            // hack to work around gcc optimisation bug
            // just expand the contents of _inc_one here:
            if(*first == '\n')
            {
               last_line = ++first;
               ++clines;
            }
            else
               ++first;
         #else         
            _inc_one(clines, last_line, first); 
         #endif
         }
         if(first == last)
            break;

         if( access::can_start(*first, _map, (unsigned char)mask_any) )
         {
            if(query_match_aux(first, last, m, e, flags, pd, &restart))
            {
               m.set_line(clines, last_line);
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
               if(m[0].second == last)
                  return cmatches;
               // update to end of what matched
               // trying to match again with match_not_null set if this
               // is a null match...
               need_init = true;
               if(first == m[0].second)
               {
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
                  if(query_match_aux(first, last, m, e, flags | match_not_null, pd, &restart))
                  {
                     m.set_line(clines, last_line);
                     ++cmatches;
                     if(foo(m) == false)
                        return cmatches;
                  }
                  else
                  {
                     need_init = false;
                     for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                         {} // dwa 10/20/2000 - warning suppression for MWCW
                     if(restart != last)
                        ++restart;
                     _skip_and_inc(clines, last_line, first, restart);
                  }
               }
               if(need_init)
               {
                  _skip_and_inc(clines, last_line, first, m[0].second);
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
               }
            }
            else
            {
               for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                   {} // dwa 10/20/2000 - warning suppression for MWCW
               if(restart != last)
                  ++restart;
               _skip_and_inc(clines, last_line, first, restart);
            }
         }
         else
            _inc_one(clines, last_line, first);
      }
   }
   break;
   case regbase::restart_line:
   {
      // do search optimised for line starts:
      while(first != last)
      {
         // find first charcter after a line break:
         --first;
         if(*first == '\n')
            --clines;
         while((first != last) && (*first != '\n'))
            ++first;
         if(first == last)
            break;
         ++first;
         if(first == last)
            break;

         ++clines;
         last_line = first;

         if( access::can_start(*first, _map, (unsigned char)mask_any) )
         {
            if(query_match_aux(first, last, m, e, flags, pd, &restart))
            {
               m.set_line(clines, last_line);
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
               if(m[0].second == last)
                  return cmatches;
               // update to end of what matched
               // trying to match again with match_not_null set if this
               // is a null match...
               need_init = true;
               if(first == m[0].second)
               {
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
                  if(query_match_aux(first, last, m, e, flags | match_not_null, pd, &restart))
                  {
                     m.set_line(clines, last_line);
                     ++cmatches;
                     if(foo(m) == false)
                        return cmatches;
                  }
                  else
                  {
                     need_init = false;
                     for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                         {} // dwa 10/20/2000 - warning suppression for MWCW
                     if(restart != last)
                        ++restart;
                     _skip_and_inc(clines, last_line, first, restart);
                  }
               }
               if(need_init)
               {
                  _skip_and_inc(clines, last_line, first, m[0].second);
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
               }
            }
            else
            {
               for(unsigned int i = 0; (restart != first) && (i < access::leading_length(e)); ++i, --restart)
                   {} // dwa 10/20/2000 - warning suppression for MWCW
               if(restart != last)
                  ++restart;
               _skip_and_inc(clines, last_line, first, restart);
            }
         }
         else
            _inc_one(clines, last_line, first);
      }
   }
   break;
   case regbase::restart_continue:
   {
      if(!leading_match)
         return cmatches;
      while(first != last)
      {
         if( access::can_start(*first, _map, (unsigned char)mask_any) )
         {
            if(query_match_aux(first, last, m, e, flags, pd, &restart))
            {
               m.set_line(clines, last_line);
               ++cmatches;
               if(foo(m) == false)
                  return cmatches;
               if(m[0].second == last)
                  return cmatches;
               // update to end of what matched
               // trying to match again with match_not_null set if this
               // is a null match...
               if(first == m[0].second)
               {
                  next_base = m[0].second;
                  pd.temp_match.init_fail(next_base, last);
                  m.init_fail(next_base, last);
                  if(query_match_aux(first, last, m, e, flags | match_not_null, pd, &restart))
                  {
                     m.set_line(clines, last_line);
                     ++cmatches;
                     if(foo(m) == false)
                        return cmatches;
                  }
                  else
                     return cmatches;  // can't continue from null match
               }
               _skip_and_inc(clines, last_line, first, m[0].second);
               next_base = m[0].second;
               pd.temp_match.init_fail(next_base, last);
               m.init_fail(next_base, last);
               continue;
            }
         }
         return cmatches;
      }
   }
   break;
   }


   // finally check trailing null string:
   if(access::first(e)->can_be_null)
   {
      if(query_match_aux(first, last, m, e, flags, pd, &restart))
      {
         m.set_line(clines, last_line);
         ++cmatches;
         if(foo(m) == false)
            return cmatches;
      }
   }

   return cmatches;
}
#if defined(BOOST_REGEX_NO_TEMPLATE_SWITCH_MERGE)
} // namespace {anon}
#endif

} // namespace re_detail

//
// proc regex_match
// returns true if the specified regular expression matches
// the whole of the input.  Fills in what matched in m.
//
template <class iterator, class Allocator, class charT, class traits, class Allocator2>
bool regex_match(iterator first, iterator last, match_results<iterator, Allocator>& m, const reg_expression<charT, traits, Allocator2>& e, unsigned flags = match_default)
{
   // prepare m for failure:
   if((flags & match_init) == 0)
   {
      m.set_size(e.mark_count(), first, last);
      m.set_base(first);
      m.set_line(1, first);
   }
   flags |= match_all; // must match all of input.
   re_detail::_priv_match_data<iterator, Allocator> pd(m, first, last, e.size());
   iterator restart;
   bool result = re_detail::query_match_aux(first, last, m, e, flags, pd, &restart);
   return result;
}
template <class iterator, class charT, class traits, class Allocator2>
bool regex_match(iterator first, iterator last, const reg_expression<charT, traits, Allocator2>& e, unsigned flags = match_default)
{
   match_results<iterator> m;
   return regex_match(first, last, m, e, flags);
}
//
// query_match convenience interfaces:
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
//
// this isn't really a partial specialisation, but template function
// overloading - if the compiler doesn't support partial specialisation
// then it really won't support this either:
template <class charT, class Allocator, class traits, class Allocator2>
inline bool regex_match(const charT* str, 
                        match_results<const charT*, Allocator>& m, 
                        const reg_expression<charT, traits, Allocator2>& e, 
                        unsigned flags = match_default)
{
   return regex_match(str, str + traits::length(str), m, e, flags);
}

template <class ST, class SA, class Allocator, class charT, class traits, class Allocator2>
inline bool regex_match(const std::basic_string<charT, ST, SA>& s, 
                 match_results<typename std::basic_string<charT, ST, SA>::const_iterator, Allocator>& m, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 unsigned flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
template <class charT, class traits, class Allocator2>
inline bool regex_match(const charT* str, 
                        const reg_expression<charT, traits, Allocator2>& e, 
                        unsigned flags = match_default)
{
   match_results<const charT*> m;
   return regex_match(str, str + traits::length(str), m, e, flags);
}

template <class ST, class SA, class charT, class traits, class Allocator2>
inline bool regex_match(const std::basic_string<charT, ST, SA>& s, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 unsigned flags = match_default)
{
   typedef typename std::basic_string<charT, ST, SA>::const_iterator iterator;
   match_results<iterator> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#else  // partial ordering
inline bool regex_match(const char* str, 
                        cmatch& m, 
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_match(str, str + regex::traits_type::length(str), m, e, flags);
}
inline bool regex_match(const char* str, 
                        const regex& e, 
                        unsigned flags = match_default)
{
   match_results<const char*> m;
   return regex_match(str, str + regex::traits_type::length(str), m, e, flags);
}
#ifndef BOOST_NO_WREGEX
inline bool regex_match(const wchar_t* str, 
                        wcmatch& m, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_match(str, str + wregex::traits_type::length(str), m, e, flags);
}
inline bool regex_match(const wchar_t* str, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   match_results<const wchar_t*> m;
   return regex_match(str, str + wregex::traits_type::length(str), m, e, flags);
}
#endif
inline bool regex_match(const std::string& s, 
                        match_results<std::string::const_iterator, regex::allocator_type>& m,
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
inline bool regex_match(const std::string& s, 
                        const regex& e, 
                        unsigned flags = match_default)
{
   match_results<std::string::const_iterator, regex::allocator_type> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#if !defined(BOOST_NO_WREGEX)
inline bool regex_match(const std::basic_string<wchar_t>& s, 
                        match_results<std::basic_string<wchar_t>::const_iterator, wregex::allocator_type>& m,
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_match(s.begin(), s.end(), m, e, flags);
}
inline bool regex_match(const std::basic_string<wchar_t>& s, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   match_results<std::basic_string<wchar_t>::const_iterator, wregex::allocator_type> m;
   return regex_match(s.begin(), s.end(), m, e, flags);
}
#endif

#endif

template <class iterator, class Allocator, class charT, class traits, class Allocator2>
bool regex_search(iterator first, iterator last, match_results<iterator, Allocator>& m, const reg_expression<charT, traits, Allocator2>& e, unsigned flags = match_default)
{
   if(e.flags() & regbase::failbit)
      return false;

   typedef typename traits::size_type traits_size_type;
   typedef typename traits::uchar_type traits_uchar_type;

   return re_detail::reg_grep2(re_detail::grep_search_predicate<iterator, Allocator>(&m), first, last, e, flags, m.allocator());
}

//
// regex_search convenience interfaces:
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
//
// this isn't really a partial specialisation, but template function
// overloading - if the compiler doesn't support partial specialisation
// then it really won't support this either:
template <class charT, class Allocator, class traits, class Allocator2>
inline bool regex_search(const charT* str, 
                        match_results<const charT*, Allocator>& m, 
                        const reg_expression<charT, traits, Allocator2>& e, 
                        unsigned flags = match_default)
{
   return regex_search(str, str + traits::length(str), m, e, flags);
}

template <class ST, class SA, class Allocator, class charT, class traits, class Allocator2>
inline bool regex_search(const std::basic_string<charT, ST, SA>& s, 
                 match_results<typename std::basic_string<charT, ST, SA>::const_iterator, Allocator>& m, 
                 const reg_expression<charT, traits, Allocator2>& e, 
                 unsigned flags = match_default)
{
   return regex_search(s.begin(), s.end(), m, e, flags);
}
#else  // partial specialisation
inline bool regex_search(const char* str, 
                        cmatch& m, 
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_search(str, str + regex::traits_type::length(str), m, e, flags);
}
#ifndef BOOST_NO_WREGEX
inline bool regex_search(const wchar_t* str, 
                        wcmatch& m, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_search(str, str + wregex::traits_type::length(str), m, e, flags);
}
#endif
inline bool regex_search(const std::string& s, 
                        match_results<std::string::const_iterator, regex::allocator_type>& m,
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_search(s.begin(), s.end(), m, e, flags);
}
#if !defined(BOOST_NO_WREGEX)
inline bool regex_search(const std::basic_string<wchar_t>& s, 
                        match_results<std::basic_string<wchar_t>::const_iterator, wregex::allocator_type>& m,
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_search(s.begin(), s.end(), m, e, flags);
}
#endif

#endif


//
// regex_grep:
// find all non-overlapping matches within the sequence first last:
//
template <class Predicate, class iterator, class charT, class traits, class Allocator>
inline unsigned int regex_grep(Predicate foo, iterator first, iterator last, const reg_expression<charT, traits, Allocator>& e, unsigned flags = match_default)
{
   return re_detail::reg_grep2(foo, first, last, e, flags, e.allocator());
}

//
// regex_grep convenience interfaces:
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
//
// this isn't really a partial specialisation, but template function
// overloading - if the compiler doesn't support partial specialisation
// then it really won't support this either:
template <class Predicate, class charT, class Allocator, class traits>
inline unsigned int regex_grep(Predicate foo, const charT* str, 
                        const reg_expression<charT, traits, Allocator>& e, 
                        unsigned flags = match_default)
{
   return regex_grep(foo, str, str + traits::length(str), e, flags);
}

template <class Predicate, class ST, class SA, class Allocator, class charT, class traits>
inline unsigned int regex_grep(Predicate foo, const std::basic_string<charT, ST, SA>& s, 
                 const reg_expression<charT, traits, Allocator>& e, 
                 unsigned flags = match_default)
{
   return regex_grep(foo, s.begin(), s.end(), e, flags);
}
#else  // partial specialisation
inline unsigned int regex_grep(bool (*foo)(const cmatch&), const char* str, 
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_grep(foo, str, str + regex::traits_type::length(str), e, flags);
}
#ifndef BOOST_NO_WREGEX
inline unsigned int regex_grep(bool (*foo)(const wcmatch&), const wchar_t* str, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_grep(foo, str, str + wregex::traits_type::length(str), e, flags);
}
#endif
inline unsigned int regex_grep(bool (*foo)(const match_results<std::string::const_iterator, regex::allocator_type>&), const std::string& s,
                        const regex& e, 
                        unsigned flags = match_default)
{
   return regex_grep(foo, s.begin(), s.end(), e, flags);
}
#if !defined(BOOST_NO_WREGEX)
inline unsigned int regex_grep(bool (*foo)(const match_results<std::basic_string<wchar_t>::const_iterator, wregex::allocator_type>&), 
                     const std::basic_string<wchar_t>& s, 
                        const wregex& e, 
                        unsigned flags = match_default)
{
   return regex_grep(foo, s.begin(), s.end(), e, flags);
}
#endif

#endif


#ifdef __BORLANDC__
  #pragma option pop
#endif

} // namespace boost

#endif   // BOOST_REGEX_MATCH_HPP


















