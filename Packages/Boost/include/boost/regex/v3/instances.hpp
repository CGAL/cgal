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
  *   FILE         instances.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Defines those template instances that are placed in the
  *                library rather than in the users object files.
  */

//
// note no include guard, we may include this multiple times:
//
#ifndef BOOST_REGEX_NO_EXTERNAL_TEMPLATES

namespace boost{

//
// this header can be included multiple times, each time with
// a different character type, BOOST_REGEX_CHAR_T must be defined
// first:
//
#ifndef BOOST_REGEX_CHAR_T
#  error "BOOST_REGEX_CHAR_T not defined"
#endif

//
// what follows is compiler specific:
//

#ifdef __BORLANDC__

#pragma option push -a8 -b -Vx -Ve -pc

#  ifndef BOOST_REGEX_INSTANTIATE
#     pragma option push -Jgx
#  endif

template class BOOST_REGEX_DECL reg_expression< BOOST_REGEX_CHAR_T >;
template class BOOST_REGEX_DECL re_detail::match_results_base<BOOST_REGEX_CHAR_T const*>;
template class BOOST_REGEX_DECL re_detail::match_results_base<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator>;
template class BOOST_REGEX_DECL match_results<BOOST_REGEX_CHAR_T const*>;
template class BOOST_REGEX_DECL match_results<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator>;


#  ifndef BOOST_REGEX_INSTANTIATE
#     pragma option pop
#  endif

#pragma option pop

#elif defined(BOOST_MSVC)

#  ifndef BOOST_REGEX_INSTANTIATE
#     define template extern template
#  endif

#pragma warning(push)
#pragma warning(disable : 4251 4231 4660)

//template class BOOST_REGEX_DECL regex_traits< BOOST_REGEX_CHAR_T >;
template class BOOST_REGEX_DECL reg_expression< BOOST_REGEX_CHAR_T >;
template class BOOST_REGEX_DECL re_detail::match_results_base<BOOST_REGEX_CHAR_T const*, BOOST_DEFAULT_ALLOCATOR(re_detail::def_alloc_param_traits<BOOST_REGEX_CHAR_T const*>::type) >;
template class BOOST_REGEX_DECL re_detail::match_results_base<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator, BOOST_DEFAULT_ALLOCATOR(re_detail::def_alloc_param_traits<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator>::type) >;
template class BOOST_REGEX_DECL match_results<BOOST_REGEX_CHAR_T const*, BOOST_DEFAULT_ALLOCATOR(re_detail::def_alloc_param_traits<BOOST_REGEX_CHAR_T const*>::type) >;
template class BOOST_REGEX_DECL match_results<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator, BOOST_DEFAULT_ALLOCATOR(re_detail::def_alloc_param_traits<std::basic_string<BOOST_REGEX_CHAR_T>::const_iterator>::type) >;


#pragma warning(pop)

#  ifdef template
#     undef template
#  endif

#elif !defined(BOOST_REGEX_HAS_DLL_RUNTIME)

//
// for each [member] function declare a full specialisation of that
// [member] function, then instantiate it in one translation unit.
// This is not guarenteed to work according to the standard, but in
// practice it should work for all compilers (unless they use a realy
// perverse name mangling convention).  Unfortunately this approach
// does *not* work for Win32 style import/export, because that can
// alter the class layout.
//

#  ifndef BOOST_REGEX_INSTANTIATE
#     define template template<>
#  endif

template reg_expression<BOOST_REGEX_CHAR_T>::reg_expression(const BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T)&);
template reg_expression<BOOST_REGEX_CHAR_T>::reg_expression(const BOOST_REGEX_CHAR_T* p, reg_expression<BOOST_REGEX_CHAR_T>::flag_type f, const BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T)& a);
template reg_expression<BOOST_REGEX_CHAR_T>::reg_expression(const BOOST_REGEX_CHAR_T* p, size_type len, reg_expression<BOOST_REGEX_CHAR_T>::flag_type f, const BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T)&);
template reg_expression<BOOST_REGEX_CHAR_T>::reg_expression(const reg_expression<BOOST_REGEX_CHAR_T>&);
template reg_expression<BOOST_REGEX_CHAR_T>::~reg_expression();
template reg_expression<BOOST_REGEX_CHAR_T>& BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::operator=(const reg_expression&);
template BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::get_allocator()const;
template bool BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::operator==(const reg_expression<BOOST_REGEX_CHAR_T>&)const;
template bool BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::operator<(const reg_expression<BOOST_REGEX_CHAR_T>&)const;
template BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::allocator()const;
template unsigned int BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::set_expression(const BOOST_REGEX_CHAR_T* p, const BOOST_REGEX_CHAR_T* end, reg_expression<BOOST_REGEX_CHAR_T>::flag_type f);
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_maps();
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_map(re_detail::re_syntax_base* node, unsigned char* _map, unsigned int* pnull, unsigned char mask, re_detail::re_syntax_base* terminal)const;
template bool BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::probe_start(re_detail::re_syntax_base* node, BOOST_REGEX_CHAR_T c, re_detail::re_syntax_base* terminal)const;
template bool BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::probe_start_null(re_detail::re_syntax_base* node, re_detail::re_syntax_base* terminal)const;
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::fixup_apply(re_detail::re_syntax_base* b, unsigned cbraces);
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::move_offsets(re_detail::re_syntax_base* j, unsigned size);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_set(const BOOST_REGEX_CHAR_T*& first, const BOOST_REGEX_CHAR_T* last);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_set_aux(re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& singles, re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T)>& ranges, re_detail::jstack<boost::uint_fast32_t, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& classes, re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& equivalents, bool isnot, const re_detail::_narrow_type&);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_set_aux(re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& singles, re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T)>& ranges, re_detail::jstack<boost::uint_fast32_t, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& classes, re_detail::jstack<reg_expression<BOOST_REGEX_CHAR_T>::traits_string_type, BOOST_DEFAULT_ALLOCATOR(BOOST_REGEX_CHAR_T) >& equivalents, bool isnot, const re_detail::_wide_type&);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::compile_set_simple(re_detail::re_syntax_base* dat, unsigned long cls, bool isnot);
template unsigned int BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::parse_inner_set(const BOOST_REGEX_CHAR_T*& first, const BOOST_REGEX_CHAR_T* last);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::add_simple(re_detail::re_syntax_base* dat, re_detail::syntax_element_type type, unsigned int size);
template re_detail::re_syntax_base* BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::add_literal(re_detail::re_syntax_base* dat, BOOST_REGEX_CHAR_T c);
template BOOST_REGEX_CHAR_T BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::parse_escape(const BOOST_REGEX_CHAR_T*& first, const BOOST_REGEX_CHAR_T* last);
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::parse_range(const BOOST_REGEX_CHAR_T*& first, const BOOST_REGEX_CHAR_T* last, unsigned& min, unsigned& max);
template bool BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::skip_space(const BOOST_REGEX_CHAR_T*& first, const BOOST_REGEX_CHAR_T* last);
template unsigned int BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::probe_restart(re_detail::re_syntax_base* dat);
template unsigned int BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::fixup_leading_rep(re_detail::re_syntax_base* dat, re_detail::re_syntax_base* end);
template void BOOST_REGEX_CALL reg_expression<BOOST_REGEX_CHAR_T>::fail(unsigned int err);

namespace re_detail{

#define iterator const BOOST_REGEX_CHAR_T*
#define Allocator match_results_base<iterator>::alloc_type
#define size_type match_results_base<iterator>::size_type

template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::set_first(iterator i);
template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::set_first(iterator i, std::size_t pos);
template match_results_base<iterator, Allocator>::match_results_base(const Allocator& a);
template Allocator BOOST_REGEX_CALL match_results_base<iterator, Allocator>::allocator()const;
template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::m_free();
template bool match_results_base<iterator, Allocator>::operator==(const match_results_base<iterator, Allocator>& that)const;
template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::set_size(size_type n);
template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::set_size(size_type n, iterator i, iterator j);
template void BOOST_REGEX_CALL match_results_base<iterator, Allocator>::cow();

#undef iterator
#undef Allocator
#undef size_type

} // namespace re_detail

#  ifdef template
#     undef template
#  endif

#endif

} // namespace boost

#endif // BOOST_REGEX_NO_EXTERNAL_TEMPLATES
 



