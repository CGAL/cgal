/* Copyright 2003-2004 Joaquín M López Muñoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/multi_index for library home page.
 */

#ifndef BOOST_MULTI_INDEX_DETAIL_ACCESS_SPECIFIER_HPP
#define BOOST_MULTI_INDEX_DETAIL_ACCESS_SPECIFIER_HPP

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

/* In those compilers that do not accept the member template friend syntax,
 * some protected and private sections might need to be specified as
 * public.
 * As per a discussion on the Boost mailing list, and pending the
 * resolution of whether BOOST_NO_MEMBER_TEMPLATE_FRIENDS should
 * apply to MSVC 8.0, I act here as if it did. The relevant
 * discussion can be found at:
 *   [boost] [config] seems like VC 8.0 needsBOOST_NO_MEMBER_TEMPLATE_FRIENDS
 *   http://lists.boost.org/MailArchives/boost/msg68369.php
 */

#if defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS) ||\
    defined(BOOST_MSVC)&&(BOOST_MSVC==1400)
#define BOOST_MULTI_INDEX_NO_MEMBER_TEMPLATE_FRIENDS
#endif

#if defined(BOOST_MULTI_INDEX_NO_MEMBER_TEMPLATE_FRIENDS)
#define BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS public
#define BOOST_MULTI_INDEX_PRIVATE_IF_MEMBER_TEMPLATE_FRIENDS public
#else
#define BOOST_MULTI_INDEX_PROTECTED_IF_MEMBER_TEMPLATE_FRIENDS protected
#define BOOST_MULTI_INDEX_PRIVATE_IF_MEMBER_TEMPLATE_FRIENDS private
#endif

/* GCC does not correctly support in-class using declarations for template
 * functions. See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=9810
 * MSVC 7.1/8.0 seem to have a similar problem, though the conditions in
 * which the error happens are not that simple. I have yet to isolate this
 * into a snippet suitable for bug reporting.
 */

#if BOOST_WORKAROUND(__GNUC__, <3)||\
    BOOST_WORKAROUND(__GNUC__,==3)&&(__GNUC_MINOR__<4)||\
    BOOST_WORKAROUND(BOOST_MSVC,==1310)||\
    BOOST_WORKAROUND(BOOST_MSVC,==1400)
#define BOOST_MULTI_INDEX_PRIVATE_IF_USING_DECL_FOR_TEMPL_FUNCTIONS public
#else
#define BOOST_MULTI_INDEX_PRIVATE_IF_USING_DECL_FOR_TEMPL_FUNCTIONS private
#endif

#endif
