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
  *   FILE         regex_compile.hpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Declares reg_expression<> member functions.  This is
  *                an internal header file, do not include directly.
  */

#ifndef BOOST_REGEX_COMPILE_HPP
#define BOOST_REGEX_COMPILE_HPP

namespace boost{
#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif
#ifdef __BORLANDC__
#pragma option push -w-8004
#endif

namespace re_detail{


template <class traits>
struct kmp_translator
{
   typedef typename traits::char_type char_type;
   bool icase;
   const traits* pt;
   kmp_translator(bool c, traits* p) : icase(c), pt(p) {}
   char_type operator()(char_type c)
   {
      return pt->translate(c, icase);
   }
};


template <class charT, class traits_type, class Allocator>
bool BOOST_REGEX_CALL re_maybe_set_member(charT c,
                                 const re_set_long* set_,
                                 const reg_expression<charT, traits_type, Allocator>& e)
{
   const charT* p = reinterpret_cast<const charT*>(set_+1);
   bool icase = e.flags() & regex_constants::icase;
   charT col = e.get_traits().translate(c, icase);
   for(unsigned int i = 0; i < set_->csingles; ++i)
   {
      if(col == *p)
         return set_->isnot ? false : true;

      while(*p)++p;
      ++p;     // skip null
   }
   return set_->isnot ? true : false;
}

} // namespace re_detail

template <class traits>
struct is_big_char
{
   typedef typename traits::uchar_type traits_uchar_type;
   typedef typename traits::size_type traits_size_type;
#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
   static bool test(char)
   { return false; }
   static bool test(unsigned char)
   { return false; }
   static bool test(signed char)
   { return false; }
#endif
   template <class charT> static bool test(charT c)
   { return (traits_size_type)(traits_uchar_type)c >= 256; }
};

template <class charT, class traits, class Allocator>
inline bool BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::can_start(charT c, const unsigned char* _map, unsigned char mask, const re_detail::_wide_type&)
{
   if(is_big_char<traits>::test(c))
      return true;
   return BOOST_REGEX_MAKE_BOOL(_map[(traits_uchar_type)c] & mask);
}

template <class charT, class traits, class Allocator>
inline bool BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::can_start(charT c, const unsigned char* _map, unsigned char mask, const re_detail::_narrow_type&)
{
   return BOOST_REGEX_MAKE_BOOL(_map[(traits_uchar_type)c] & mask);
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::reg_expression(const Allocator& a)
    : regbase(), data(a), pkmp(0), error_code_(REG_EMPTY), _expression(0)
{
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::reg_expression(const charT* p, flag_type f, const Allocator& a)
    : data(a), pkmp(0), error_code_(REG_EMPTY), _expression(0)
{
   set_expression(p, f | regex_constants::use_except);
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::reg_expression(const charT* p1, const charT* p2, flag_type f, const Allocator& a)
    : data(a), pkmp(0), error_code_(REG_EMPTY), _expression(0)
{
    set_expression(p1, p2, f | regex_constants::use_except);
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::reg_expression(const charT* p, size_type len, flag_type f, const Allocator& a)
    : data(a), pkmp(0), error_code_(REG_EMPTY), _expression(0)
{
    set_expression(p, p + len, f | regex_constants::use_except);
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::reg_expression(const reg_expression<charT, traits, Allocator>& e)
   : regbase(e), data(e.allocator()), pkmp(0), error_code_(REG_EMPTY), _expression(0)
{
   //
   // we do a deep copy only if e is a valid expression, otherwise fail.
   //
   if(e.error_code() == 0)
   {
      const charT* pe = e.expression();
      set_expression(pe, pe + e._expression_len, e.flags() | regex_constants::use_except);
   }
   else
   {
      _flags = e.flags() & ~(regex_constants::use_except);
      fail(e.error_code());
   }
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>::~reg_expression()
{
   if(pkmp)
      re_detail::kmp_free(pkmp, data.allocator());
}

template <class charT, class traits, class Allocator>
reg_expression<charT, traits, Allocator>& BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::operator=(const reg_expression<charT, traits, Allocator>& e)
{
   //
   // we do a deep copy only if e is a valid expression, otherwise fail.
   //
   if(this == &e) return *this;
   _flags = use_except;
   fail(e.error_code());
   if(error_code() == 0)
      set_expression(e._expression, e._expression + e._expression_len, e.flags() | regex_constants::use_except);
   return *this;
}

template <class charT, class traits, class Allocator>
int BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compare(const reg_expression<charT, traits, Allocator>& e)const
{
   if(_flags != e.flags())
      return _flags - e.flags();
   return str().compare(e.str());
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::swap(reg_expression& that)throw()
{
   traits_inst.swap(that.traits_inst);
   data.swap(that.data);
   static_cast<regbase&>(*this).swap(that);

   std::swap(_restart_type, that._restart_type);
   std::swap(marks, that.marks);
   std::swap(repeats, that.repeats);
   std::swap(startmap, that.startmap);
   std::swap(_expression_len, that._expression_len);
   std::swap(_leading_len, that._leading_len);
   std::swap(_leading_string, that._leading_string);
   std::swap(_leading_string_len, that._leading_string_len);
   std::swap(pkmp, that.pkmp);
   std::swap(error_code_, that.error_code_);
   std::swap(_expression, that._expression);
}

template <class charT, class traits, class Allocator>
Allocator BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::allocator()const
{
    return data.allocator();
}

template <class charT, class traits, class Allocator>
Allocator BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::get_allocator()const
{
    return data.allocator();
}

template <class charT, class traits, class Allocator>
unsigned int BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::parse_inner_set(const charT*& arg_first, const charT* arg_last)
{
   //
   // we have an inner [...] construct
   //
   jm_assert(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) == traits_type::syntax_open_set);
   const charT* base = arg_first;
   while( (arg_first != arg_last)
      && (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) != traits_type::syntax_close_set) )
         ++arg_first;
   if(arg_first == arg_last)
      return 0;
   ++arg_first;
   if((arg_first-base) < 5)
      return 0;
   if(*(base+1) != *(arg_first-2))
      return 0;
   unsigned int result = traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*(base+1));
   if((result == traits_type::syntax_colon) && ((arg_first-base) == 5))
   {
      unsigned type = traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*(base+2));
      if((type == traits_type::syntax_left_word) || (type == traits_type::syntax_right_word))
         return type;
   }
   return ((result == traits_type::syntax_colon) || (result == traits_type::syntax_dot) || (result == traits_type::syntax_equal)) ? result : 0;
}


template <class charT, class traits, class Allocator>
bool BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::skip_space(const charT*& arg_first, const charT* arg_last)
{
   //
   // returns true if we get to arg_last:
   //
   while((arg_first != arg_last) && (traits_inst.is_class(*arg_first, traits_type::char_class_space) == true))
   {
      ++arg_first;
   }
   return arg_first == arg_last;
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::parse_range(const charT*& ptr, const charT* arg_end, unsigned& min, unsigned& max)
{
   //
   // we have {x} or {x,} or {x,y} NB no spaces inside braces
   // anything else is illegal
   // On input ptr points to "{"
   //
   ++ptr;
   if(skip_space(ptr, arg_end))
   {
      fail(REG_EBRACE);
      return;
   }
   if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) != traits_type::syntax_digit)
   {
      fail(REG_BADBR);
      return;
   }
   min = traits_inst.toi(ptr, arg_end, 10);
   if(skip_space(ptr, arg_end))
   {
      fail(REG_EBRACE);
      return;
   }
   if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) == traits_type::syntax_comma)
   {
      //we have a second interval:
      ++ptr;
      if(skip_space(ptr, arg_end))
      {
         fail(REG_EBRACE);
         return;
      }
      if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) == traits_type::syntax_digit)
         max = traits_inst.toi(ptr, arg_end, 10);
      else
         max = (unsigned)-1;
   }
   else
      max = min;

   // validate input:
   if(skip_space(ptr, arg_end))
   {
      fail(REG_EBRACE);
      return;
   }
   if(max < min)
   {
      fail(REG_ERANGE);
      return;
   }
   if(_flags & bk_braces)
   {
      if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) != traits_type::syntax_slash)
      {
         fail(REG_BADBR);
         return;
      }
      else
      {
         // back\ is OK now check the }
         ++ptr;
         if((ptr == arg_end) || (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) != traits_type::syntax_close_brace))
         {
            fail(REG_BADBR);
            return;
         }
      }
   }
   else if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) != traits_type::syntax_close_brace)
   {
      fail(REG_BADBR);
      return;
   }
}

template <class charT, class traits, class Allocator>
charT BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::parse_escape(const charT*& arg_first, const charT* arg_last)
{
   charT c(*arg_first);
   traits_size_type c_unsigned = (traits_size_type)(traits_uchar_type)*arg_first;
   // this is only used for the switch(), but cannot be folded in
   // due to a bug in Comeau 4.2.44beta3
   traits_size_type syntax = traits_inst.syntax_type(c_unsigned);
   switch(syntax)
   {
   case traits_type::syntax_a:
      c = '\a';
      ++arg_first;
      break;
   case traits_type::syntax_f:
      c = '\f';
      ++arg_first;
      break;
   case traits_type::syntax_n:
      c = '\n';
      ++arg_first;
      break;
   case traits_type::syntax_r:
      c = '\r';
      ++arg_first;
      break;
   case traits_type::syntax_t:
      c = '\t';
      ++arg_first;
      break;
   case traits_type::syntax_v:
      c = '\v';
      ++arg_first;
      break;
   case traits_type::syntax_x:
      ++arg_first;
      if(arg_first == arg_last)
      {
         fail(REG_EESCAPE);
         break;
      }
      // maybe have \x{ddd}
      if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*arg_first)) == traits_type::syntax_open_brace)
      {
         ++arg_first;
         if(arg_first == arg_last)
         {
            fail(REG_EESCAPE);
            break;
         }
         if(traits_inst.is_class(*arg_first, traits_type::char_class_xdigit) == false)
         {
            fail(REG_BADBR);
            break;
         }
         c = (charT)traits_inst.toi(arg_first, arg_last, -16);
         if((arg_first == arg_last) || (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)(*arg_first)) != traits_type::syntax_close_brace))
         {
            fail(REG_BADBR);
         }
         ++arg_first;
         break;
      }
      else
      {
         if(traits_inst.is_class(*arg_first, traits_type::char_class_xdigit) == false)
         {
            fail(REG_BADBR);
            break;
         }
         c = (charT)traits_inst.toi(arg_first, arg_last, -16);
      }
      break;
   case traits_type::syntax_c:
      ++arg_first;
      if(arg_first == arg_last)
      {
         fail(REG_EESCAPE);
         break;
      }
      if(((traits_uchar_type)(*arg_first) < (traits_uchar_type)'@')
            || ((traits_uchar_type)(*arg_first) > (traits_uchar_type)127) )
      {
         fail(REG_EESCAPE);
         return (charT)0;
      }
      c = (charT)((traits_uchar_type)(*arg_first) - (traits_uchar_type)'@');
      ++arg_first;
      break;
   case traits_type::syntax_e:
      c = (charT)27;
      ++arg_first;
      break;
   case traits_type::syntax_digit:
      c = (charT)traits_inst.toi(arg_first, arg_last, -8);
      break;
   default:
      //c = *arg_first;
      ++arg_first;
   }
   return c;
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_maps()
{
   re_detail::re_syntax_base* record = static_cast<re_detail::re_syntax_base*>(data.data());
   // always compile the first _map:
   std::memset(startmap, 0, 256);
   record->can_be_null = 0;
   compile_map(record, startmap, 0, re_detail::mask_all);

   while(record->type != re_detail::syntax_element_match)
   {
      if((record->type == re_detail::syntax_element_alt) || (record->type == re_detail::syntax_element_rep))
      {
         std::memset(&(static_cast<re_detail::re_jump*>(record)->_map), 0, 256);
         record->can_be_null = 0;
         compile_map(record->next.p, static_cast<re_detail::re_jump*>(record)->_map, &(record->can_be_null), re_detail::mask_take, static_cast<re_detail::re_jump*>(record)->alt.p);
         compile_map(static_cast<re_detail::re_jump*>(record)->alt.p, static_cast<re_detail::re_jump*>(record)->_map, &(record->can_be_null), re_detail::mask_skip);
         if(record->type == re_detail::syntax_element_rep)
         {
            re_detail::re_repeat* rep = static_cast<re_detail::re_repeat*>(record);
            // set whether this is a singleton repeat or not:
            if(rep->next.p->next.p->next.p == rep->alt.p)
            {
               switch(rep->next.p->type)
               {
               case re_detail::syntax_element_wild:
                  rep->type = re_detail::syntax_element_dot_rep;
                  break;
               case re_detail::syntax_element_literal:
                  rep->type = re_detail::syntax_element_char_rep;
                  break;
               case re_detail::syntax_element_set:
                  rep->type = re_detail::syntax_element_short_set_rep;
                  break;
               case re_detail::syntax_element_long_set:
                  if(static_cast<re_detail::re_set_long*>(rep->next.p)->singleton)
                     rep->type = re_detail::syntax_element_long_set_rep;
                  break;
               default:
                  break;
               }
            }
         }
      }
      else
      {
         record->can_be_null = 0;
         compile_map(record, 0, &(record->can_be_null), re_detail::mask_all);
      }
      record = record->next.p;
   }
   record->can_be_null = re_detail::mask_all;
}

template <class charT, class traits, class Allocator>
bool BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::probe_start(
                        re_detail::re_syntax_base* node, charT cc, re_detail::re_syntax_base* terminal) const
{
   unsigned int c;

   switch(node->type)
   {
   case re_detail::syntax_element_startmark:
      if(static_cast<const re_detail::re_brace*>(node)->index == -1)
      {
         return probe_start(node->next.p->next.p, cc, terminal)
            && probe_start(static_cast<const re_detail::re_jump*>(node->next.p)->alt.p, cc, terminal);
      }
      else if(static_cast<const re_detail::re_brace*>(node)->index == -3)
      {
         return probe_start(node->next.p->next.p, cc, terminal);
      }
      // doesn't tell us anything about the next character, so:
      return probe_start(node->next.p, cc, terminal);
   case re_detail::syntax_element_endmark:
      if(static_cast<const re_detail::re_brace*>(node)->index == -3)
      {
         return true;
      }
      // fall through:
   case re_detail::syntax_element_start_line:
   case re_detail::syntax_element_word_boundary:
   case re_detail::syntax_element_buffer_start:
   case re_detail::syntax_element_restart_continue:
      // doesn't tell us anything about the next character, so:
      return probe_start(node->next.p, cc, terminal);
   case re_detail::syntax_element_literal:
      // only the first character of the literal can match:
      // note these have already been translated:
      if(*reinterpret_cast<charT*>(static_cast<re_detail::re_literal*>(node)+1) == traits_inst.translate(cc, (_flags & regex_constants::icase)))
         return true;
      return false;
   case re_detail::syntax_element_end_line:
      // next character (if there is one!) must be a newline:
      if(traits_inst.is_separator(traits_inst.translate(cc, (_flags & regex_constants::icase))))
         return true;
      return false;
   case re_detail::syntax_element_wild:
      return true;
   case re_detail::syntax_element_match:
      return true;
   case re_detail::syntax_element_within_word:
   case re_detail::syntax_element_word_start:
      return traits_inst.is_class(traits_inst.translate(cc, (_flags & regex_constants::icase)), traits_type::char_class_word);
   case re_detail::syntax_element_word_end:
      // what follows must not be a word character,
      return traits_inst.is_class(traits_inst.translate(cc, (_flags & regex_constants::icase)), traits_type::char_class_word) ? false : true;
   case re_detail::syntax_element_buffer_end:
      // we can be null, nothing must follow,
      // NB we assume that this is followed by
      // re_detail::syntax_element_match, if its not then we can
      // never match anything anyway!!
      return false;
   case re_detail::syntax_element_soft_buffer_end:
      // we can be null, only newlines must follow,
      // NB we assume that this is followed by
      // re_detail::syntax_element_match, if its not then we can
      // never match anything anyway!!
      return traits_inst.is_separator(traits_inst.translate(cc, (_flags & regex_constants::icase)));
   case re_detail::syntax_element_backref:
      // there's no easy way to determine this
      // which is not to say it can't be done!
      // for now:
      return true;
   case re_detail::syntax_element_long_set:
      // we can not be null,
      // we need to add already translated values in the set
      // to values in the _map
      return re_detail::re_maybe_set_member(cc, static_cast<const re_detail::re_set_long*>(node), *this) || (re_detail::re_is_set_member(static_cast<const charT*>(&cc), static_cast<const charT*>(&cc+1), static_cast<re_detail::re_set_long*>(node), *this) != &cc);
   case re_detail::syntax_element_set:
      // set all the elements that are set in corresponding set:
      c = (traits_size_type)(traits_uchar_type)traits_inst.translate(cc, (_flags & regex_constants::icase));
      return static_cast<re_detail::re_set*>(node)->_map[c] != 0;
   case re_detail::syntax_element_jump:
      if(static_cast<re_detail::re_jump*>(node)->alt.p < node)
      {
         // backwards jump,
         // caused only by end of repeat section, we'll treat this
         // the same as a match, because the sub-expression has matched.
         if(node->next.p == terminal)
            return true; // null repeat - we can always take this
         else
         {
            //
            // take the jump, add in fix for the fact that if the
            // repeat that we're jumping to has non-zero minimum count
            // then we need to add in the possiblity that we could still
            // skip that repeat.
            re_detail::re_syntax_base* next = static_cast<re_detail::re_jump*>(node)->alt.p;
            bool b = probe_start(next, cc, terminal);
            if((next->type == re_detail::syntax_element_rep) && (static_cast<re_detail::re_repeat*>(next)->min != 0))
            {
               b = b || probe_start(static_cast<re_detail::re_jump*>(next)->alt.p, cc, terminal);
            }
            return b;
         }
      }
      else
         // take the jump and compile:
         return probe_start(static_cast<re_detail::re_jump*>(node)->alt.p, cc, terminal);
   case re_detail::syntax_element_alt:
      // we need to take the OR of the two alternatives:
      return probe_start(static_cast<re_detail::re_jump*>(node)->alt.p, cc, terminal) || probe_start(node->next.p, cc, terminal);
   case re_detail::syntax_element_rep:
   case re_detail::syntax_element_char_rep:
   case re_detail::syntax_element_dot_rep:
   case re_detail::syntax_element_long_set_rep:
   case re_detail::syntax_element_short_set_rep:
      // we need to take the OR of the two alternatives
      if(static_cast<re_detail::re_repeat*>(node)->min == 0)
         return probe_start(node->next.p, cc, static_cast<re_detail::re_jump*>(node)->alt.p) || probe_start(static_cast<re_detail::re_jump*>(node)->alt.p, cc, terminal);
      else
         return probe_start(node->next.p, cc, static_cast<re_detail::re_jump*>(node)->alt.p);
   case re_detail::syntax_element_combining:
      return !traits_inst.is_combining(traits_inst.translate(cc, (_flags & regex_constants::icase)));
   }
   return false;
}

template <class charT, class traits, class Allocator>
bool BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::probe_start_null(re_detail::re_syntax_base* node, re_detail::re_syntax_base* terminal)const
{
   switch(node->type)
   {
   case re_detail::syntax_element_endmark:
      if(static_cast<const re_detail::re_brace*>(node)->index == -3)
      {
         return true;
      }
      // fall through:
   case re_detail::syntax_element_startmark:
   case re_detail::syntax_element_start_line:
   case re_detail::syntax_element_word_boundary:
   case re_detail::syntax_element_buffer_start:
   case re_detail::syntax_element_restart_continue:
   case re_detail::syntax_element_end_line:
   case re_detail::syntax_element_word_end:
      // doesn't tell us anything about the next character, so:
      return probe_start_null(node->next.p, terminal);
   case re_detail::syntax_element_match:
   case re_detail::syntax_element_buffer_end:
   case re_detail::syntax_element_soft_buffer_end:
   case re_detail::syntax_element_backref:
      return true;
   case re_detail::syntax_element_jump:
      if(static_cast<re_detail::re_jump*>(node)->alt.p < node)
      {
         // backwards jump,
         // caused only by end of repeat section, we'll treat this
         // the same as a match, because the sub-expression has matched.
         // this is only caused by NULL repeats as in "(a*)*" or "(\<)*"
         // these are really nonsensence and make the matching code much
         // harder, it would be nice to get rid of them altogether.
         if(node->next.p == terminal)
            return true;
         else
            return probe_start_null(static_cast<re_detail::re_jump*>(node)->alt.p, terminal);
      }
      else
         // take the jump and compile:
         return probe_start_null(static_cast<re_detail::re_jump*>(node)->alt.p, terminal);
   case re_detail::syntax_element_alt:
      // we need to take the OR of the two alternatives:
      return probe_start_null(static_cast<re_detail::re_jump*>(node)->alt.p, terminal) || probe_start_null(node->next.p, terminal);
   case re_detail::syntax_element_rep:
      // only need to consider skipping the repeat:
      return probe_start_null(static_cast<re_detail::re_jump*>(node)->alt.p, terminal);
   default:
      break;
   }
   return false;
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_map(
                        re_detail::re_syntax_base* node, unsigned char* _map,
                        unsigned int* pnull, unsigned char mask, re_detail::re_syntax_base* terminal)const
{
   if(_map)
   {
      for(unsigned int i = 0; i < 256; ++i)
      {
         if(probe_start(node, (charT)i, terminal))
            _map[i] |= mask;
      }
   }
   if(pnull && probe_start_null(node, terminal))
      *pnull |= mask;
}
  
template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::move_offsets(re_detail::re_syntax_base* j, unsigned arg_size)
{
# ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable: 4127)
#endif
   // move all offsets starting with j->link forward by arg_size
   // called after an insert:
   j = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(data.data()) + j->next.i);
   while(true)
   {
      switch(j->type)
      {
      case re_detail::syntax_element_rep:
         static_cast<re_detail::re_jump*>(j)->alt.i += arg_size;
         j->next.i += arg_size;
         break;
      case re_detail::syntax_element_jump:
      case re_detail::syntax_element_alt:
         static_cast<re_detail::re_jump*>(j)->alt.i += arg_size;
         j->next.i += arg_size;
         break;
      default:
         j->next.i += arg_size;
         break;
      }
      if(j->next.i == arg_size)
         break;
      j = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(data.data()) + j->next.i);
   }
# ifdef BOOST_MSVC
#  pragma warning(pop)
#endif
}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_set_simple(re_detail::re_syntax_base* dat, unsigned long cls, bool isnot)
{
   typedef typename re_detail::is_byte<charT>::width_type width_type;
   re_detail::jstack<traits_string_type, Allocator> singles(64, data.allocator());
   re_detail::jstack<traits_string_type, Allocator> ranges(64, data.allocator());
   re_detail::jstack<boost::uint_fast32_t, Allocator> classes(64, data.allocator());
   re_detail::jstack<traits_string_type, Allocator> equivalents(64, data.allocator());
   if(_flags & regbase::icase)
   {
      if((cls == traits_type::char_class_upper) || (cls == traits_type::char_class_lower))
      {
         cls = traits_type::char_class_alpha;
      }
   }
   classes.push(cls);
   if(dat)
   {
      data.align();
      dat->next.i = data.size();
   }
   return compile_set_aux(singles, ranges, classes, equivalents, isnot, width_type());
}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_set(const charT*& arg_first, const charT* arg_last)
{
   re_detail::jstack<traits_string_type, Allocator> singles(64, data.allocator());
   re_detail::jstack<traits_string_type, Allocator> ranges(64, data.allocator());
   re_detail::jstack<boost::uint_fast32_t, Allocator> classes(64, data.allocator());
   re_detail::jstack<traits_string_type, Allocator> equivalents(64, data.allocator());
   bool has_digraphs = false;
   jm_assert(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) == traits_type::syntax_open_set);
   ++arg_first;
   bool started = false;
   bool done = false;
   bool isnot = false;

   enum last_type
   {
      last_single,
      last_none,
      last_dash
   };

   unsigned l = last_none;
   traits_string_type s;

   while((arg_first != arg_last) && !done)
   {
      traits_size_type c = (traits_size_type)(traits_uchar_type)*arg_first;
      // this is only used for the switch(), but cannot be folded in
      // due to a bug in Comeau 4.2.44beta3
      traits_size_type syntax = traits_inst.syntax_type(c);
      switch(syntax)
      {
      case traits_type::syntax_caret:
         if(!started && !isnot)
         {
            isnot = true;
         }
         else
         {
            s = (charT)c;
            goto char_set_literal;
         }
         break;
      case traits_type::syntax_open_set:
      {
         if((_flags & char_classes) == 0)
         {
            s = (charT)c;
            goto char_set_literal;
         }
         // check to see if we really have a class:
         const charT* base = arg_first;
         // this is only used for the switch(), but cannot be folded in
         // due to a bug in Comeau 4.2.44beta3
    unsigned int inner_set = parse_inner_set(arg_first, arg_last);
         switch(inner_set)
         {
         case traits_type::syntax_colon:
            {
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               boost::uint_fast32_t id = traits_inst.lookup_classname(base+2, arg_first-2);
               if(_flags & regex_constants::icase)
               {
                  if((id == traits_type::char_class_upper) || (id == traits_type::char_class_lower))
                  {
                     id = traits_type::char_class_alpha;
                  }
               }
               if(id == 0)
               {
                  fail(REG_ECTYPE);
                  return 0;
               }
               classes.push(id);
               started = true;
               l = last_none;
            }
            break;
         case traits_type::syntax_dot:
            //
            // we have a collating element [.collating-name.]
            //
            if(traits_inst.lookup_collatename(s, base+2, arg_first-2))
            {
               --arg_first;
               if(s.size() > 1)
                  has_digraphs = true;
               if(s.size())goto char_set_literal;
            }
            fail(REG_ECOLLATE);
            return 0;
         case traits_type::syntax_equal:
            //
            // we have an equivalence class [=collating-name=]
            //
            if(traits_inst.lookup_collatename(s, base+2, arg_first-2))
            {
               std::size_t len = s.size();
               if(len)
               {
                  unsigned i = 0;
                  while(i < len)
                  {
                     s[i] = traits_inst.translate(s[i], (_flags & regex_constants::icase));
                     ++i;
                  }
                  traits_string_type s2;
                  traits_inst.transform_primary(s2, s);
                  equivalents.push(s2);
                  started = true;
                  l = last_none;
                  break;
               }
            }
            fail(REG_ECOLLATE);
            return 0;
         case traits_type::syntax_left_word:
            if((started == false) && (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) == traits_type::syntax_close_set))
            {
               ++arg_first;
               return add_simple(0, re_detail::syntax_element_word_start);
            }
            fail(REG_EBRACK);
            return 0;
         case traits_type::syntax_right_word:
            if((started == false) && (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) == traits_type::syntax_close_set))
            {
               ++arg_first;
               return add_simple(0, re_detail::syntax_element_word_end);
            }
            fail(REG_EBRACK);
            return 0;
         default:
            if(started == false)
            {
               unsigned int t = traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*(base+1));
               if((t != traits_type::syntax_colon) && (t != traits_type::syntax_dot) && (t != traits_type::syntax_equal))
               {
                  arg_first = base;
                  s = (charT)c;
                  goto char_set_literal;
               }
            }
            fail(REG_EBRACK);
            return 0;
         }
         if(arg_first == arg_last)
         {
            fail(REG_EBRACK);
            return 0;
         }
         continue;
      }
      case traits_type::syntax_close_set:
         if(started == false)
         {
            s = (charT)c;
            goto char_set_literal;
         }
         done = true;
         break;
      case traits_type::syntax_dash:
         if(!started)
         {
            s = (charT)c;
            goto char_set_literal;
         }
         ++arg_first;
         if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*arg_first) == traits_type::syntax_close_set)
         {
            --arg_first;
            s = (charT)c;
            goto char_set_literal;
         }
         if((singles.empty() == true) || (l != last_single))
         {
            fail(REG_ERANGE);
            return 0;
         }
         ranges.push(singles.peek());
         if(singles.peek().size() <= 1)  // leave digraphs and ligatures in place
            singles.pop();
         l = last_dash;
         continue;
      case traits_type::syntax_slash:
         if(_flags & regex_constants::escape_in_lists)
         {
            ++arg_first;
            if(arg_first == arg_last)
               continue;
            /*traits_size_type*/ c = (traits_size_type)(traits_uchar_type)*arg_first;
            // this is only used for the switch(), but cannot be folded in
            // due to a bug in Comeau 4.2.44beta3
            traits_size_type syntax4 = traits_inst.syntax_type(c);
            switch(syntax4)
            {
            case traits_type::syntax_w:
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               classes.push(traits_type::char_class_word);
               started = true;
               l = last_none;
               ++arg_first;
               continue;
            case traits_type::syntax_d:
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               classes.push(traits_type::char_class_digit);
               started = true;
               l = last_none;
               ++arg_first;
               continue;
            case traits_type::syntax_s:
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               classes.push(traits_type::char_class_space);
               started = true;
               l = last_none;
               ++arg_first;
               continue;
            case traits_type::syntax_l:
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               classes.push(traits_type::char_class_lower);
               started = true;
               l = last_none;
               ++arg_first;
               continue;
            case traits_type::syntax_u:
               if(l == last_dash)
               {
                  fail(REG_ERANGE);
                  return 0;
               }
               classes.push(traits_type::char_class_upper);
               started = true;
               l = last_none;
               ++arg_first;
               continue;
            case traits_type::syntax_W:
            case traits_type::syntax_D:
            case traits_type::syntax_S:
            case traits_type::syntax_U:
            case traits_type::syntax_L:
               fail(REG_EESCAPE);
               return 0;
            default:
               c = parse_escape(arg_first, arg_last);
               --arg_first;
               s = (charT)c;
               goto char_set_literal;
            }
         }
         else
         {
            s = (charT)c;
            goto char_set_literal;
         }
      default:
         s = (charT)c;
         char_set_literal:
         unsigned i = 0;
         // get string length to stop us going past the end of string (DWA)
         std::size_t len = s.size();
         while(i < len)
         {
            s[i] = traits_inst.translate(s[i], (_flags & regex_constants::icase));
            ++i;
         }
         started = true;
         if(l == last_dash)
         {
            ranges.push(s);
            l = last_none;
            if(s.size() > 1)   // add ligatures to singles list as well
               singles.push(s);
         }
         else
         {
            singles.push(s);
            l = last_single;
         }
      }
      ++arg_first;
   }
   if(!done)
      return 0;

   typedef typename re_detail::is_byte<charT>::width_type width_type;
   
   re_detail::re_syntax_base* result;
   if(has_digraphs)
      result = compile_set_aux(singles, ranges, classes, equivalents, isnot, re_detail::_wide_type());
   else
      result = compile_set_aux(singles, ranges, classes, equivalents, isnot, width_type());
   #ifdef __BORLANDC__
   // delayed throw:
   if((result == 0) && (_flags & regex_constants::use_except))
      fail(error_code());
   #endif
   return result;
}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_set_aux(re_detail::jstack<traits_string_type, Allocator>& singles, re_detail::jstack<traits_string_type, Allocator>& ranges, re_detail::jstack<boost::uint_fast32_t, Allocator>& classes, re_detail::jstack<traits_string_type, Allocator>& equivalents, bool isnot, const re_detail::_wide_type&)
{
   size_type base = data.size();
   data.extend(sizeof(re_detail::re_set_long));
   unsigned int csingles = 0;
   unsigned int cranges = 0;
   boost::uint_fast32_t cclasses = 0;
   unsigned int cequivalents = 0;
   bool nocollate_state = !(flags() & regex_constants::collate);
   bool singleton = true;

   while(singles.empty() == false)
   {
      ++csingles;
      const traits_string_type& s = singles.peek();
      std::size_t len = (s.size() + 1) * sizeof(charT);
      if(len > sizeof(charT) * 2)
         singleton = false;
      std::memcpy(reinterpret_cast<charT*>(data.extend(len)), s.c_str(), len);
      singles.pop();
   }
   while(ranges.empty() == false)
   {
      traits_string_type c1, c2;
      if(nocollate_state)
         c1 = ranges.peek();
      else
         traits_inst.transform(c1, ranges.peek());
      ranges.pop();
      if(nocollate_state)
         c2 = ranges.peek();
      else
         traits_inst.transform(c2, ranges.peek());
      ranges.pop();
      if(c1 < c2)
      {
         // for some reason bc5 crashes when throwing exceptions
         // from here - probably an EH-compiler bug, but hard to
         // be sure...
         // delay throw to later:
         #ifdef __BORLANDC__
         boost::uint_fast32_t f = _flags;
         _flags &= ~regex_constants::use_except;
         #endif
         fail(REG_ERANGE);
         #ifdef __BORLANDC__
         _flags = f;
         #endif
         return 0;
      }
      ++cranges;
      std::size_t len = (re_detail::re_strlen(c1.c_str()) + 1) * sizeof(charT);
      std::memcpy(data.extend(len), c1.c_str(), len);
      len = (re_detail::re_strlen(c2.c_str()) + 1) * sizeof(charT);
      std::memcpy(data.extend(len), c2.c_str(), len);
   }
   while(classes.empty() == false)
   {
      cclasses |= classes.peek();
      classes.pop();
   }
   while(equivalents.empty() == false)
   {
      ++cequivalents;
      const traits_string_type& s = equivalents.peek();
      std::size_t len = (re_detail::re_strlen(s.c_str()) + 1) * sizeof(charT);
      std::memcpy(reinterpret_cast<charT*>(data.extend(len)), s.c_str(), len);
      equivalents.pop();
   }

   re_detail::re_set_long* dat = reinterpret_cast<re_detail::re_set_long*>(reinterpret_cast<unsigned char*>(data.data()) + base);
   dat->type = re_detail::syntax_element_long_set;
   dat->csingles = csingles;
   dat->cranges = cranges;
   dat->cclasses = cclasses;
   dat->cequivalents = cequivalents;
   dat->isnot = isnot;
   dat->next.i = 0;
   dat->singleton = isnot ? true : singleton;
   return dat;
}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::compile_set_aux(re_detail::jstack<traits_string_type, Allocator>& singles, re_detail::jstack<traits_string_type, Allocator>& ranges, re_detail::jstack<boost::uint_fast32_t, Allocator>& classes, re_detail::jstack<traits_string_type, Allocator>& equivalents, bool isnot, const re_detail::_narrow_type&)
{
   re_detail::re_set* dat = reinterpret_cast<re_detail::re_set*>(data.extend(sizeof(re_detail::re_set)));
   std::memset(dat, 0, sizeof(re_detail::re_set));

   while(singles.empty() == false)
   {
      dat->_map[(traits_size_type)(traits_uchar_type)*(singles.peek().c_str())] = re_detail::mask_all;
      singles.pop();
   }
   while(ranges.empty() == false)
   {
      traits_string_type c1, c2, c3, c4;

      if((flags() & regex_constants::collate) == 0)
         c1 = ranges.peek();
      else
         traits_inst.transform(c1, ranges.peek());
      ranges.pop();
      if((flags() & regex_constants::collate) == 0)
         c2 = ranges.peek();
      else
         traits_inst.transform(c2, ranges.peek());
      ranges.pop();

      if(c1 < c2)
      {
         // for some reason bc5 crashes when throwing exceptions
         // from here - probably an EH-compiler bug, but hard to
         // be sure...
         // delay throw to later:
         #ifdef __BORLANDC__
         boost::uint_fast32_t f = _flags;
         _flags &= ~regex_constants::use_except;
         #endif
         fail(REG_ERANGE);
         #ifdef __BORLANDC__
         _flags = f;
         #endif
         return 0;
      }
      for(unsigned int i = 0; i < 256; ++i)
      {
         c4 = (charT)i;
         if((flags() & regex_constants::collate) == 0)
            c3 = c4;
         else
            traits_inst.transform(c3, c4);
         if((c3 <= c1) && (c3 >= c2))
            dat->_map[i] = re_detail::mask_all;
      }
   }
   while(equivalents.empty() == false)
   {
      traits_string_type c1, c2;
      for(unsigned int i = 0; i < 256; ++i)
      {
         c2 = (charT)i;
         traits_inst.transform_primary(c1, c2);
         if(c1 == equivalents.peek())
            dat->_map[i] = re_detail::mask_all;
      }
      equivalents.pop();
   }

   boost::uint_fast32_t l_flags = 0;
   while(classes.empty() == false)
   {
      l_flags |= classes.peek();
      classes.pop();
   }
   if(l_flags)
   {
      for(unsigned int i = 0; i < 256; ++i)
      {
         if(traits_inst.is_class(charT(i), l_flags))
            dat->_map[(traits_uchar_type)traits_inst.translate((charT)i, (_flags & regex_constants::icase))] = re_detail::mask_all;
      }
   }

   if(isnot)
   {
      for(unsigned int i = 0; i < 256; ++i)
      {
         dat->_map[i] = !dat->_map[i];
      }
   }

   dat->type = re_detail::syntax_element_set;
   dat->next.i = 0;
   return dat;
}

#ifndef __CODEGUARD__
// this must not be inline when Borland's codeguard support is turned
// on, otherwise we _will_ get surious codeguard errors...
inline
#endif
 re_detail::re_syntax_base* add_offset(void* base, std::ptrdiff_t off)
{
   return reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(base) + off);
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::fixup_apply(re_detail::re_syntax_base* b, unsigned cbraces)
{
   typedef typename boost::detail::rebind_allocator<bool, Allocator>::type b_alloc;
   
   register unsigned char* base = reinterpret_cast<unsigned char*>(b);
   register re_detail::re_syntax_base* ptr = b;
   bool* pb = 0;
   b_alloc a(data.allocator());
#ifndef BOOST_NO_EXCEPTIONS
   try
   {
#endif
      pb = a.allocate(cbraces);
      BOOST_REGEX_NOEH_ASSERT(pb)
      for(unsigned i = 0; i < cbraces; ++i)
         pb[i] = false;

      repeats = 0;

      while(ptr->next.i)
      {
         switch(ptr->type)
         {
         case re_detail::syntax_element_rep:
            jm_assert(data.size() > static_cast<re_detail::re_jump*>(ptr)->alt.i);
            static_cast<re_detail::re_jump*>(ptr)->alt.p = add_offset(base, static_cast<re_detail::re_jump*>(ptr)->alt.i);
#ifdef BOOST_REGEX_DEBUG
            if((re_detail::padding_mask & reinterpret_cast<int>(static_cast<re_detail::re_jump*>(ptr)->alt.p)) && (static_cast<re_detail::re_jump*>(ptr)->alt.p != b))
            {
               jm_trace("padding mis-aligment in repeat jump to object type: " << static_cast<re_detail::re_jump*>(ptr)->alt.p->type)
               //jm_assert(0 == (padding_mask & (int)((re_detail::re_jump*)ptr)->alt.p));
            }
#endif
            static_cast<re_detail::re_repeat*>(ptr)->id = repeats;
            ++repeats;
            goto rebase;
         case re_detail::syntax_element_jump:
         case re_detail::syntax_element_alt:
            jm_assert(data.size() > static_cast<re_detail::re_jump*>(ptr)->alt.i);
            static_cast<re_detail::re_jump*>(ptr)->alt.p = add_offset(base, static_cast<re_detail::re_jump*>(ptr)->alt.i);
#ifdef BOOST_REGEX_DEBUG
            if((re_detail::padding_mask & reinterpret_cast<int>(static_cast<re_detail::re_jump*>(ptr)->alt.p) && (static_cast<re_detail::re_jump*>(ptr)->alt.p != b)))
            {
               jm_trace("padding mis-aligment in alternation jump to object type: " << static_cast<re_detail::re_jump*>(ptr)->alt.p->type)
               //jm_assert(0 == (padding_mask & (int)((re_detail::re_jump*)ptr)->alt.p));
            }
#endif
            goto rebase;
         case re_detail::syntax_element_backref:
            if((static_cast<re_detail::re_brace*>(ptr)->index >= (int)cbraces) || (pb[static_cast<re_detail::re_brace*>(ptr)->index] == false) )
            {
               fail(REG_ESUBREG);
               a.deallocate(pb, cbraces);
               return;
            }
            goto rebase;
         case re_detail::syntax_element_endmark:
            if(static_cast<re_detail::re_brace*>(ptr)->index > 0)
               pb[static_cast<re_detail::re_brace*>(ptr)->index] = true;
            goto rebase;
         default:
            rebase:
            jm_assert(data.size() > ptr->next.i);
            ptr->next.p = add_offset(base, ptr->next.i);
#ifdef BOOST_REGEX_DEBUG
            if((re_detail::padding_mask & (int)(ptr->next.p)) && (static_cast<re_detail::re_jump*>(ptr)->alt.p != b))
            {
               jm_trace("padding mis-alignment in next record of type " << ptr->next.p->type)
               jm_assert(0 == (re_detail::padding_mask & (int)(ptr->next.p)));
            }
#endif
            ptr = ptr->next.p;
         }
      }
      a.deallocate(pb, cbraces);
      pb = 0;
#ifndef BOOST_NO_EXCEPTIONS
   }
   catch(...)
   {
      if(pb)
         a.deallocate(pb, cbraces);
      throw;
   }
#endif
}


template <class charT, class traits, class Allocator>
unsigned int BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::set_expression(const charT* arg_first, const charT* arg_last, flag_type f)
{
# ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable: 4127)
#endif
#ifdef __OpenBSD__ 
   // strxfrm not working on OpenBSD?? 
   f &= ~regex_constants::collate; 
#endif 

   if(arg_first == expression())
   {
      traits_string_type s(arg_first, arg_last);
      return set_expression(s.c_str(), s.c_str() + s.size(), f);
   }
   typedef typename traits_type::sentry sentry_t;
   sentry_t sent(traits_inst);
   if(sent){

   const charT* base = arg_first;
   data.clear();
   _flags = f;
   fail(REG_NOERROR);  // clear any error
   _leading_len = 0; // set this to non-zero if there are any backrefs, we'll refer to it later...

   if(arg_first >= arg_last)
   {
      fail(REG_EMPTY);
      return error_code();
   }

   const charT* ptr = arg_first;
   marks = 0;
   re_detail::jstack<std::size_t, Allocator> mark(64, data.allocator());
   re_detail::jstack<int, Allocator> markid(64, data.allocator());
   std::size_t last_mark_popped = 0;
   register traits_size_type c;
   register re_detail::re_syntax_base* dat;

   unsigned rep_min = 0;
   unsigned rep_max = 0;

   //
   // set up header:
   //
   ++marks;
   dat = 0;

   if(_flags & regex_constants::literal)
   {
      while(ptr != arg_last)
      {
         dat = add_literal(dat, traits_inst.translate(*ptr, (_flags & regex_constants::icase)));
         ++ptr;
      }
   }

   while (ptr < arg_last)
   {
      c = (traits_size_type)(traits_uchar_type)*ptr;
      // this is only used for the switch(), but cannot be folded in
      // due to a bug in Comeau 4.2.44beta3
      traits_size_type syntax = traits_inst.syntax_type(c);
      switch(syntax)
      {
      case traits_type::syntax_open_bracket:
         if(_flags & bk_parens)
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         open_bracked_jump:
         // extend:
         dat = add_simple(dat, re_detail::syntax_element_startmark, sizeof(re_detail::re_brace));
         if(_flags & nosubs)
         {
            markid.push(0);
            static_cast<re_detail::re_brace*>(dat)->index = 0;
         }
         else
         {
            markid.push(marks);
            static_cast<re_detail::re_brace*>(dat)->index = marks++;
         }
         mark.push(data.index(dat));
         ++ptr;
         //
         // check for perl like (?...) extention syntax
         c = (traits_size_type)(traits_uchar_type)*ptr;
         if(((_flags & (bk_parens|perlex)) == perlex) && (traits_type::syntax_question == traits_inst.syntax_type(c)))
         {
            ++ptr;
            c = (traits_size_type)(traits_uchar_type)*ptr;
            // this is only used for the switch(), but cannot be folded in
            // due to a bug in Comeau 4.2.44beta3
            traits_size_type syntax2 = traits_inst.syntax_type(c);
            switch(syntax2)
            {
            case traits_type::syntax_colon:
               static_cast<re_detail::re_brace*>(dat)->index = 0;
               if((_flags & nosubs) == 0)
                  --marks;
               markid.pop();
               markid.push(0);
               ++ptr;
               continue;
            case traits_type::syntax_equal:
               static_cast<re_detail::re_brace*>(dat)->index = -1;
               markid.pop();
               markid.push(-1);
               common_forward_assert:
               if((_flags & nosubs) == 0)
                  --marks;
               ++ptr;
               // extend:
               dat = add_simple(dat, re_detail::syntax_element_jump, re_detail::re_jump_size);
               data.align();
               //
               // we don't know what value to put here yet,
               // use an arbitrarily large value for now
               // and check it later:
               static_cast<re_detail::re_jump*>(dat)->alt.i = INT_MAX/2;
               mark.push(data.size() - re_detail::re_jump_size);
               continue;
            case traits_type::syntax_right_word:
               static_cast<re_detail::re_brace*>(dat)->index = -3;
               markid.pop();
               markid.push(-3);
               goto common_forward_assert;
            case traits_type::syntax_not:
               static_cast<re_detail::re_brace*>(dat)->index = -2;
               markid.pop();
               markid.push(-2);
               goto common_forward_assert;
            case traits_type::syntax_hash:
               // comment just skip it:
               static_cast<re_detail::re_brace*>(dat)->index = 0;
               if((_flags & nosubs) == 0)
                  --marks;
               markid.pop();
               mark.pop();
               do{
                  ++ptr;
                  c = (traits_size_type)(traits_uchar_type)*ptr;
               }while(traits_type::syntax_close_bracket != traits_inst.syntax_type(c));
               ++ptr;
               continue;
            default:
               //
               // error, return to standard parsing and let that handle the error:
               --ptr;
               continue;
            }
         }
         break;
      case traits_type::syntax_close_bracket:
         if(_flags & bk_parens)
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         
         close_bracked_jump:
         if(dat)
         {
            data.align();
            dat->next.i = data.size();
         }

         if(mark.empty())
         {
            fail(REG_EPAREN);
            return error_code();
         }
         // see if we have an empty alternative:
         if(mark.peek() == data.index(dat) )
         {
            re_detail::re_syntax_base* para = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(data.data()) + mark.peek());
            if(para->type == re_detail::syntax_element_jump)
            {
               fail(REG_EMPTY);
               return error_code();
            }
         }

         // pop any pushed alternatives and set the target arg_last destination:
         dat = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<unsigned char*>(data.data()) + mark.peek());
         while(dat->type == re_detail::syntax_element_jump)
         {
            static_cast<re_detail::re_jump*>(dat)->alt.i = data.size();
            mark.pop();
            if(mark.empty())
            {
               fail(REG_EPAREN);
               return error_code();
            }
            dat = reinterpret_cast<re_detail::re_jump*>(reinterpret_cast<unsigned char*>(data.data()) + mark.peek());
         }

         dat = add_simple(0, re_detail::syntax_element_endmark, sizeof(re_detail::re_brace));
         static_cast<re_detail::re_brace*>(dat)->index = markid.peek();
         markid.pop();
         last_mark_popped = mark.peek();
         mark.pop();
         ++ptr;
         break;
      case traits_type::syntax_char:
         dat = add_literal(dat, (charT)c);
         ++ptr;
         break;
      case traits_type::syntax_slash:
      {
         if(++ptr == arg_last)
         {
            fail(REG_EESCAPE);
            return error_code();
         }
         c = (traits_size_type)(traits_uchar_type)*ptr;
         // this is only used for the switch(), but cannot be folded in
         // due to a bug in Comeau 4.2.44beta3
         traits_size_type syntax3 = traits_inst.syntax_type(c);
         switch(syntax3)
         {
         case traits_type::syntax_open_bracket:
            if(_flags & bk_parens)
               goto open_bracked_jump;
            break;
         case traits_type::syntax_close_bracket:
            if(_flags & bk_parens)
               goto close_bracked_jump;
            break;
         case traits_type::syntax_plus:
            if((_flags & bk_plus_qm) && ((_flags & limited_ops) == 0))
            {
               rep_min = 1;
               rep_max = (unsigned)-1;
               goto repeat_jump;
            }
            break;
         case traits_type::syntax_question:
            if((_flags & bk_plus_qm) && ((_flags & limited_ops) == 0))
            {
               rep_min = 0;
               rep_max = 1;
               goto repeat_jump;
            }
            break;
         case traits_type::syntax_or:
            if(((_flags & bk_vbar) == 0) || (_flags & limited_ops))
               break;
            goto alt_string_jump;
         case traits_type::syntax_open_brace:
            if( ((_flags & bk_braces) == 0) || ((_flags & intervals) == 0))
               break;

            // we have {x} or {x,} or {x,y}:
            parse_range(ptr, arg_last, rep_min, rep_max);
            goto repeat_jump;

         case traits_type::syntax_digit:
            if(_flags & bk_refs)
            {
               // update previous:
               int i = traits_inst.toi((charT)c);
               if(i == 0)
               {
                  // we can have \025 which means take char whose
                  // code is 25 (octal), so parse string:
                  c = traits_inst.toi(ptr, arg_last, -8);
                  --ptr;
                  break;
               }
               dat = add_simple(dat, re_detail::syntax_element_backref, sizeof(re_detail::re_brace));
               static_cast<re_detail::re_brace*>(dat)->index = i;
               ++ptr;
               _leading_len = 1;
               continue;
            }
            break;
         case traits_type::syntax_b:     // re_detail::syntax_element_word_boundary
            dat = add_simple(dat, re_detail::syntax_element_word_boundary);
            ++ptr;
            continue;
         case traits_type::syntax_B:
            dat = add_simple(dat, re_detail::syntax_element_within_word);
            ++ptr;
            continue;
         case traits_type::syntax_left_word:
            dat = add_simple(dat, re_detail::syntax_element_word_start);
            ++ptr;
            continue;
         case traits_type::syntax_right_word:
            dat = add_simple(dat, re_detail::syntax_element_word_end);
            ++ptr;
            continue;
         case traits_type::syntax_w:     //re_detail::syntax_element_word_char
            dat = compile_set_simple(dat, traits_type::char_class_word);
            ++ptr;
            continue;
         case traits_type::syntax_W:
            dat = compile_set_simple(dat, traits_type::char_class_word, true);
            ++ptr;
            continue;
         case traits_type::syntax_d:     //re_detail::syntax_element_word_char
            dat = compile_set_simple(dat, traits_type::char_class_digit);
            ++ptr;
            continue;
         case traits_type::syntax_D:
            dat = compile_set_simple(dat, traits_type::char_class_digit, true);
            ++ptr;
            continue;
         case traits_type::syntax_s:     //re_detail::syntax_element_word_char
            dat = compile_set_simple(dat, traits_type::char_class_space);
            ++ptr;
            continue;
         case traits_type::syntax_S:
            dat = compile_set_simple(dat, traits_type::char_class_space, true);
            ++ptr;
            continue;
         case traits_type::syntax_l:     //re_detail::syntax_element_word_char
            dat = compile_set_simple(dat, traits_type::char_class_lower);
            ++ptr;
            continue;
         case traits_type::syntax_L:
            dat = compile_set_simple(dat, traits_type::char_class_lower, true);
            ++ptr;
            continue;
         case traits_type::syntax_u:     //re_detail::syntax_element_word_char
            dat = compile_set_simple(dat, traits_type::char_class_upper);
            ++ptr;
            continue;
         case traits_type::syntax_U:
            dat = compile_set_simple(dat, traits_type::char_class_upper, true);
            ++ptr;
            continue;
         case traits_type::syntax_Q:
            ++ptr;
            while(true)
            {
               if(ptr == arg_last)
               {
                  fail(REG_EESCAPE);
                  return error_code();
               }
               if(traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) == traits_type::syntax_slash)
               {
                  ++ptr;
                  if((ptr != arg_last) && (traits_inst.syntax_type((traits_size_type)(traits_uchar_type)*ptr) == traits_type::syntax_E))
                     break;
                  else
                  {
                     dat = add_literal(dat, *(ptr-1));
                     continue;
                  }
               }
               dat = add_literal(dat, *ptr);
               ++ptr;
            }
            ++ptr;
            continue;
         case traits_type::syntax_C:
            dat = add_simple(dat, re_detail::syntax_element_wild);
            ++ptr;
            continue;
         case traits_type::syntax_X:
            dat = add_simple(dat, re_detail::syntax_element_combining);
            ++ptr;
            continue;
         case traits_type::syntax_Z:
            dat = add_simple(dat, re_detail::syntax_element_soft_buffer_end);
            ++ptr;
            continue;
         case traits_type::syntax_G:
            dat = add_simple(dat, re_detail::syntax_element_restart_continue);
            ++ptr;
            continue;
         case traits_type::syntax_start_buffer:
            dat = add_simple(dat, re_detail::syntax_element_buffer_start);
            ++ptr;
            continue;
         case traits_type::syntax_end_buffer:
            dat = add_simple(dat, re_detail::syntax_element_buffer_end);
            ++ptr;
            continue;
         default:
            c = (traits_size_type)(traits_uchar_type)parse_escape(ptr, arg_last);
            dat = add_literal(dat, (charT)c);
            continue;
         }
         dat = add_literal(dat, (charT)c);
         ++ptr;
         break;
      }
      case traits_type::syntax_dollar:
         dat = add_simple(dat, re_detail::syntax_element_end_line, sizeof(re_detail::re_syntax_base));
         ++ptr;
         continue;
      case traits_type::syntax_caret:
         dat = add_simple(dat, re_detail::syntax_element_start_line, sizeof(re_detail::re_syntax_base));
         ++ptr;
         continue;
      case traits_type::syntax_dot:
         dat = add_simple(dat, re_detail::syntax_element_wild, sizeof(re_detail::re_syntax_base));
         ++ptr;
         continue;
      case traits_type::syntax_star:
         rep_min = 0;
         rep_max = (unsigned)-1;

         repeat_jump:
         {
          std::ptrdiff_t offset;
            if(dat == 0)
            {
               fail(REG_BADRPT);
               return error_code();
            }
            switch(dat->type)
            {
            case re_detail::syntax_element_endmark:
               offset = last_mark_popped;
               break;
            case re_detail::syntax_element_literal:
               if(static_cast<re_detail::re_literal*>(dat)->length > 1)
               {
                  // update previous:
                  charT lit = *reinterpret_cast<charT*>(reinterpret_cast<char*>(dat) + sizeof(re_detail::re_literal) + ((static_cast<re_detail::re_literal*>(dat)->length-1)*sizeof(charT)));
                  --static_cast<re_detail::re_literal*>(dat)->length;
                  dat = add_simple(dat, re_detail::syntax_element_literal, sizeof(re_detail::re_literal) + sizeof(charT));
                  static_cast<re_detail::re_literal*>(dat)->length = 1;
                  *reinterpret_cast<charT*>(static_cast<re_detail::re_literal*>(dat)+1) = lit;
               }
               offset = reinterpret_cast<char*>(dat) - reinterpret_cast<char*>(data.data());
               break;
            case re_detail::syntax_element_backref:
            case re_detail::syntax_element_long_set:
            case re_detail::syntax_element_set:
            case re_detail::syntax_element_wild:
            case re_detail::syntax_element_combining:
               // we're repeating a single item:
               offset = reinterpret_cast<char*>(dat) - reinterpret_cast<char*>(data.data());
               break;
            default:
               fail(REG_BADRPT);
               return error_code();
            }
            data.align();
            dat->next.i = data.size();
            //unsigned pos = (char*)dat - (char*)data.data();

            // add the trailing jump:
            dat = add_simple(dat, re_detail::syntax_element_jump, re_detail::re_jump_size);
            static_cast<re_detail::re_jump*>(dat)->alt.i = 0;

            // now insert the leading repeater:
            dat = static_cast<re_detail::re_syntax_base*>(data.insert(offset, re_detail::re_repeater_size));
            dat->next.i = (reinterpret_cast<char*>(dat) - reinterpret_cast<char*>(data.data())) + re_detail::re_repeater_size;
            dat->type = re_detail::syntax_element_rep;
            static_cast<re_detail::re_repeat*>(dat)->alt.i = data.size();
            static_cast<re_detail::re_repeat*>(dat)->min = rep_min;
            static_cast<re_detail::re_repeat*>(dat)->max = rep_max;
            static_cast<re_detail::re_repeat*>(dat)->leading = false;
            static_cast<re_detail::re_repeat*>(dat)->greedy = true;
            move_offsets(dat, re_detail::re_repeater_size);
            ++ptr;
            //
            // now check to see if we have a non-greedy repeat:
            if((ptr != arg_last) && (_flags & (perlex | limited_ops | bk_plus_qm | bk_braces)) == perlex)
            {
               c = (traits_size_type)(traits_uchar_type)*ptr;
               if(traits_type::syntax_question == traits_inst.syntax_type(c))
               {
                  // OK repeat is non-greedy:
                  static_cast<re_detail::re_repeat*>(dat)->greedy = false;
                  ++ptr;
               }
            }
            dat = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(data.data()) + data.size() - re_detail::re_jump_size);
            static_cast<re_detail::re_repeat*>(dat)->alt.i = offset;
            continue;
         }
      case traits_type::syntax_plus:
         if(_flags & (bk_plus_qm | limited_ops))
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         rep_min = 1;
         rep_max = (unsigned)-1;
         goto repeat_jump;
      case traits_type::syntax_question:
         if(_flags & (bk_plus_qm | limited_ops))
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         rep_min = 0;
         rep_max = 1;
         goto repeat_jump;
      case traits_type::syntax_open_set:
         // update previous:
         if(dat)
         {
            data.align();
            dat->next.i = data.size();
         }
         // extend:
         dat = compile_set(ptr, arg_last);
         if(dat == 0)
         {
            if((_flags & regex_constants::failbit) == 0)
               fail(REG_EBRACK);
            return error_code();
         }
         break;
      case traits_type::syntax_or:
      {
         if(_flags & (bk_vbar | limited_ops))
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }

         alt_string_jump:

         // update previous:
         if(dat == 0)
         {
            // start of pattern can't have empty "|"
            fail(REG_EMPTY);
            return error_code();
         }
         // see if we have an empty alternative:
         if(mark.empty() == false)
            if(mark.peek() == data.index(dat))
            {
               fail(REG_EMPTY);
               return error_code();
            }
         // extend:
         dat = add_simple(dat, re_detail::syntax_element_jump, re_detail::re_jump_size);
         data.align();
         //
         // we don't know what value to put here yet,
         // use an arbitrarily large value for now
         // and check it later (TODO!)
         static_cast<re_detail::re_jump*>(dat)->alt.i = INT_MAX/2;

         // now work out where to insert:
         std::size_t offset = 0;
         if(mark.empty() == false)
         {
            // we have a '(' or '|' to go back to:
            offset = mark.peek();
            re_detail::re_syntax_base* base2 = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<unsigned char*>(data.data()) + offset);
            offset = base2->next.i;
         }
         re_detail::re_jump* j = static_cast<re_detail::re_jump*>(data.insert(offset, re_detail::re_jump_size));
         j->type = re_detail::syntax_element_alt;
         j->next.i = offset + re_detail::re_jump_size;
         j->alt.i = data.size();
         move_offsets(j, re_detail::re_jump_size);
         dat = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<unsigned char*>(data.data()) + data.size() - re_detail::re_jump_size);
         mark.push(data.size() - re_detail::re_jump_size);
         ++ptr;
         break;
      }
      case traits_type::syntax_open_brace:
         if((_flags & bk_braces) || ((_flags & intervals) == 0))
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         // we have {x} or {x,} or {x,y}:
         parse_range(ptr, arg_last, rep_min, rep_max);
         goto repeat_jump;
      case traits_type::syntax_newline:
         if(_flags & newline_alt)
            goto alt_string_jump;
         dat = add_literal(dat, (charT)c);
         ++ptr;
         continue;
      case traits_type::syntax_close_brace:
         if(_flags & bk_braces)
         {
            dat = add_literal(dat, (charT)c);
            ++ptr;
            continue;
         }
         fail(REG_BADPAT);
         return error_code();
      default:
         dat = add_literal(dat, (charT)c);
         ++ptr;
         break;
      }  // switch
   }     // while

   //
   // update previous:
   if(dat)
   {
      data.align();
      dat->next.i = data.size();
   }

   // see if we have an empty alternative:
   if(mark.empty() == false)
      if(mark.peek() == data.index(dat) )
      {
         re_detail::re_syntax_base* para = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<char*>(data.data()) + mark.peek());
         if(para->type == re_detail::syntax_element_jump)
         {
            fail(REG_EMPTY);
            return error_code();
         }
      }
   //
   // set up tail:
   //
   if(mark.empty() == false)
   {
      // pop any pushed alternatives and set the target arg_last destination:
      dat = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<unsigned char*>(data.data()) + mark.peek());
      while(dat->type == re_detail::syntax_element_jump)
      {
         static_cast<re_detail::re_jump*>(dat)->alt.i = data.size();
         mark.pop();
         if(mark.empty() == true)
            break;
         dat = reinterpret_cast<re_detail::re_jump*>(reinterpret_cast<unsigned char*>(data.data()) + mark.peek());
      }
   }

   dat = static_cast<re_detail::re_brace*>(data.extend(sizeof(re_detail::re_syntax_base)));
   dat->type = re_detail::syntax_element_match;
   dat->next.i = 0;

   if(mark.empty() == false)
   {
      fail(REG_EPAREN);
      return error_code();
   }

   //
   // allocate space for start _map:
   startmap = reinterpret_cast<unsigned char*>(data.extend(256 + ((arg_last - base + 1) * sizeof(charT))));
   //
   // and copy the expression we just compiled:
   _expression = reinterpret_cast<charT*>(reinterpret_cast<char*>(startmap) + 256);
   _expression_len = arg_last - base;
   std::memcpy(_expression, base, _expression_len * sizeof(charT));
   *(_expression + _expression_len) = charT(0);

   //
   // now we need to apply fixups to the array
   // so that we can use pointers and not indexes
   fixup_apply(static_cast<re_detail::re_syntax_base*>(data.data()), marks);

   // check for error during fixup:
   if(_flags & regex_constants::failbit)
      return error_code();

   //
   // finally compile the maps so that we can make intelligent choices
   // whenever we encounter an alternative:
   compile_maps();
   if(pkmp)
   {
      re_detail::kmp_free(pkmp, data.allocator());
      pkmp = 0;
   }
   re_detail::re_syntax_base* sbase = static_cast<re_detail::re_syntax_base*>(data.data());
   _restart_type = probe_restart(sbase);
   _leading_len = fixup_leading_rep(sbase, 0);
   if((sbase->type == re_detail::syntax_element_literal) && (sbase->next.p->type == re_detail::syntax_element_match))
   {
      _restart_type = restart_fixed_lit;
      if(0 == pkmp)
      {
         charT* p1 = reinterpret_cast<charT*>(reinterpret_cast<char*>(sbase) + sizeof(re_detail::re_literal));
         charT* p2 = p1 + static_cast<re_detail::re_literal*>(sbase)->length;
         pkmp = re_detail::kmp_compile(p1, p2, charT(), re_detail::kmp_translator<traits>(_flags&regex_constants::icase, &traits_inst), data.allocator());
      }
   }
   return error_code();

   } // sentry
   return REG_EMPTY;

# ifdef BOOST_MSVC
#  pragma warning(pop)
#endif

}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::add_simple(re_detail::re_syntax_base* dat, re_detail::syntax_element_type type, unsigned int arg_size)
{
   if(dat)
   {
      data.align();
      dat->next.i = data.size();
   }
   if(arg_size < sizeof(re_detail::re_syntax_base))
      arg_size = sizeof(re_detail::re_syntax_base);
   dat = static_cast<re_detail::re_syntax_base*>(data.extend(arg_size));
   dat->type = type;
   dat->next.i = 0;
   return dat;
}

template <class charT, class traits, class Allocator>
re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::add_literal(re_detail::re_syntax_base* dat, charT c)
{
   if(dat && (dat->type == re_detail::syntax_element_literal))
   {
      // add another charT to the list:
      std::ptrdiff_t pos = reinterpret_cast<unsigned char*>(dat) - reinterpret_cast<unsigned char*>(data.data());
      *reinterpret_cast<charT*>(data.extend(sizeof(charT))) = traits_inst.translate(c, (_flags & regex_constants::icase));
      dat = reinterpret_cast<re_detail::re_syntax_base*>(reinterpret_cast<unsigned char*>(data.data()) + pos);
      ++(static_cast<re_detail::re_literal*>(dat)->length);
   }
   else
   {
      // extend:
      dat = add_simple(dat, re_detail::syntax_element_literal, sizeof(re_detail::re_literal) + sizeof(charT));
      static_cast<re_detail::re_literal*>(dat)->length = 1;
      *reinterpret_cast<charT*>(reinterpret_cast<re_detail::re_literal*>(dat)+1) = traits_inst.translate(c, (_flags & regex_constants::icase));
   }
   return dat;
}

template <class charT, class traits, class Allocator>
unsigned int BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::probe_restart(re_detail::re_syntax_base* dat)
{
   switch(dat->type)
   {
   case re_detail::syntax_element_startmark:
   case re_detail::syntax_element_endmark:
      if(static_cast<const re_detail::re_brace*>(dat)->index == -2) 
         return regbase::restart_any; 
      return probe_restart(dat->next.p);
   case re_detail::syntax_element_start_line:
      return regbase::restart_line;
   case re_detail::syntax_element_word_start:
      return regbase::restart_word;
   case re_detail::syntax_element_buffer_start:
      return regbase::restart_buf;
   case re_detail::syntax_element_restart_continue:
      return regbase::restart_continue;
   default:
      return regbase::restart_any;
   }
}

template <class charT, class traits, class Allocator>
unsigned int BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::fixup_leading_rep(re_detail::re_syntax_base* dat, re_detail::re_syntax_base* arg_end)
{
   unsigned int len = 0;
   if((_restart_type >= restart_word) || (_restart_type <= restart_continue))
      return 0;
   bool leading_lit = arg_end ? false : true;
   while(dat != arg_end)
   {
      switch(dat->type)
      {
      case re_detail::syntax_element_literal:
         len += static_cast<re_detail::re_literal*>(dat)->length;
         if((leading_lit) && (static_cast<re_detail::re_literal*>(dat)->length > 2))
         {
            // we can do a literal search for the leading literal string
            // using Knuth-Morris-Pratt (or whatever), and only then check for 
            // matches.  We need a decent length string though to make it
            // worth while.
            _leading_string = reinterpret_cast<charT*>(reinterpret_cast<char*>(dat) + sizeof(re_detail::re_literal));
            _leading_string_len = static_cast<re_detail::re_literal*>(dat)->length;
            _restart_type = restart_lit;
            leading_lit = false;
            const charT* p1 = _leading_string;
            const charT* p2 = _leading_string + _leading_string_len;
            pkmp = re_detail::kmp_compile(p1, p2, charT(), re_detail::kmp_translator<traits>(_flags&regex_constants::icase, &traits_inst), data.allocator());
         }
         leading_lit = false;
         break;
      case re_detail::syntax_element_wild:
         ++len;
         leading_lit = false;
         break;
      case re_detail::syntax_element_match:
         return len;
      case re_detail::syntax_element_backref:
      //case re_detail::syntax_element_jump:
      case re_detail::syntax_element_alt:
      case re_detail::syntax_element_combining:
         return 0;
      case re_detail::syntax_element_long_set:
      {
         // we need to verify that there are no multi-character
         // collating elements inside the repeat:
         if(!static_cast<re_detail::re_set_long*>(dat)->singleton)
            return 0;
         ++len;
         leading_lit = false;
         break;
      }
      case re_detail::syntax_element_set:
         ++len;
         leading_lit = false;
         break;
      case re_detail::syntax_element_rep:
      case re_detail::syntax_element_dot_rep:
      case re_detail::syntax_element_char_rep:
      case re_detail::syntax_element_short_set_rep: 
      case re_detail::syntax_element_long_set_rep:
         if((len == 0) && (_leading_len == 0) && (1 == fixup_leading_rep(dat->next.p, static_cast<re_detail::re_repeat*>(dat)->alt.p) ))
         {
            static_cast<re_detail::re_repeat*>(dat)->leading = leading_lit;
            return len;
         }
         return len;
     case re_detail::syntax_element_startmark: 
         if(static_cast<const re_detail::re_brace*>(dat)->index == -2) 
            return 0; 
         // fall through: 
      default:
         break;
      }
      dat = dat->next.p;
   }
   return len;
}

template <class charT, class traits, class Allocator>
void BOOST_REGEX_CALL reg_expression<charT, traits, Allocator>::fail(unsigned int err)
{
   error_code_  = err;
   if(err)
   {
      _flags |= regex_constants::failbit;
#ifndef BOOST_NO_EXCEPTIONS
      if(_flags & regex_constants::use_except)
      {
         re_detail::raise_error(traits_inst, err);
      }
#endif
   }
   else
      _flags &= ~regex_constants::failbit;
}

#ifdef __BORLANDC__
#pragma option pop
#endif
#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

} // namespace boost


#endif   // BOOST_REGEX_COMPILE_HPP











