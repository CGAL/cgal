// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Eric Berberich    <ericb@post.tau.ac.il>

#ifndef CGAL_ARR_TAGS_H
#define CGAL_ARR_TAGS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <boost/type_traits.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/has_xxx.hpp>

/*! \file
 * Definition of the tags for the arrangement package.
 */

namespace CGAL {

struct Arr_boundary_side_tag {};
struct Arr_oblivious_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_not_oblivious_side_tag : public virtual Arr_boundary_side_tag {};
struct Arr_open_side_tag : public virtual Arr_not_oblivious_side_tag {};
struct Arr_closed_side_tag : public virtual Arr_not_oblivious_side_tag {};
struct Arr_contracted_side_tag : public virtual Arr_not_oblivious_side_tag {};
struct Arr_identified_side_tag : public virtual Arr_not_oblivious_side_tag {};

BOOST_MPL_HAS_XXX_TRAIT_DEF(Left_side_category)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Bottom_side_category)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Top_side_category)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Right_side_category)

namespace internal {

//! type to provide left side tag (is oblivious if not existing)
template <typename Traits_, bool B>
struct Get_left_side_category {};

template <typename Traits_>
struct Get_left_side_category<Traits_, true> {
  typedef typename Traits_::Left_side_category          Category;
};

template <typename Traits_>
struct Get_left_side_category<Traits_, false> {
  typedef Arr_oblivious_side_tag                        Category;
};

template <typename Traits_>
struct Arr_complete_left_side_category {
public:
  typedef Traits_                                       Traits;
  typedef typename
  Get_left_side_category<Traits,
                         has_Left_side_category<Traits>::value>::Category
    Category;
};

template <typename GeometryTraits_2, bool b>
struct Validate_left_side_category {};

template <typename GeometryTraits_2>
struct Validate_left_side_category<GeometryTraits_2, true> {
  template <typename T>
  void missing__Left_side_category() {}
};

template <typename GeometryTraits_2>
struct Validate_left_side_category<GeometryTraits_2, false> {
  template <typename T>
  void missing__Left_side_category()
  { T missing__Left_side_category__assuming__Arr_oblivious_side_tag__instead; }
};


//! type to provide bottom side tag (is oblivious if not existing)
template <typename Traits_, bool B>
struct Get_bottom_side_category {};

template <typename Traits_>
struct Get_bottom_side_category<Traits_, true> {
  typedef typename Traits_::Bottom_side_category        Category;
};

template <typename Traits_>
struct Get_bottom_side_category<Traits_, false> {
  typedef Arr_oblivious_side_tag                        Category;
};

template <typename Traits_>
struct Arr_complete_bottom_side_category {
public:
  typedef Traits_                                       Traits;
  typedef typename
  Get_bottom_side_category<Traits,
                           has_Bottom_side_category<Traits>::value>::Category
    Category;
};

template <typename GeometryTraits_2, bool b>
struct Validate_bottom_side_category {};

template <typename GeometryTraits_2>
struct Validate_bottom_side_category<GeometryTraits_2, true> {
  template <typename T>
  void missing__Bottom_side_category() {}
};

template <typename GeometryTraits_2>
struct Validate_bottom_side_category<GeometryTraits_2, false> {
  template <typename T>
  void missing__Bottom_side_category()
  { T missing__Bottom_side_category__assuming__Arr_oblivious_side_tag__instead; }
};

//! type to provide top side tag (is oblivious if not existing)
template <typename Traits_, bool B>
struct Get_top_side_category {
};

template <typename Traits_>
struct Get_top_side_category<Traits_, true> {
  typedef typename Traits_::Top_side_category           Category;
};

template <typename Traits_>
struct Get_top_side_category<Traits_, false> {
  typedef Arr_oblivious_side_tag                        Category;
};

template <typename Traits_>
struct Arr_complete_top_side_category {
public:
  typedef Traits_                                       Traits;
  typedef typename
  Get_top_side_category<Traits, has_Top_side_category<Traits>::value>::Category
    Category;
};

template <typename GeometryTraits_2, bool b>
struct Validate_top_side_category {};

template <typename GeometryTraits_2>
struct Validate_top_side_category<GeometryTraits_2, true> {
  template <typename T>
  void missing__Top_side_category() {}
};

template <typename GeometryTraits_2>
struct Validate_top_side_category<GeometryTraits_2, false> {
  template <typename T>
  void missing__Top_side_category()
  { T missing__Top_side_category__assuming__Arr_oblivious_side_tag__instead; }
};

//! type to provide right side tag (is oblivious if not existing)
template <typename Traits_, bool B>
struct Get_right_side_category {};

template <typename Traits_>
struct Get_right_side_category<Traits_, true> {
  typedef typename Traits_::Right_side_category         Category;
};

template <typename Traits_>
struct Get_right_side_category<Traits_, false> {
  typedef Arr_oblivious_side_tag                        Category;
};

template <typename Traits_>
struct Arr_complete_right_side_category {
public:
  typedef Traits_                                       Traits;

  typedef typename
  Get_right_side_category<Traits,
                          has_Right_side_category<Traits>::value>::Category
    Category;
};

template <typename GeometryTraits_2, bool b>
struct Validate_right_side_category {};

template <typename GeometryTraits_2>
struct Validate_right_side_category<GeometryTraits_2, true> {
  template <typename T>
  void missing__Right_side_category() {}
};

template <typename GeometryTraits_2>
struct Validate_right_side_category<GeometryTraits_2, false> {
  template <typename T>
  void missing__Right_side_category()
  { T missing__Right_side_category__assuming__Arr_oblivious_side_tag__instead; }
};

} // namespace internal

struct Arr_boundary_cond_tag{};
struct Arr_all_sides_oblivious_tag : public virtual Arr_boundary_cond_tag{};
struct Arr_not_all_sides_oblivious_tag : public virtual Arr_boundary_cond_tag{};

struct Arr_has_identified_side_tag :
    public virtual Arr_not_all_sides_oblivious_tag{};
struct Arr_has_contracted_side_tag :
    public virtual Arr_not_all_sides_oblivious_tag{};
struct Arr_has_closed_side_tag :
    public virtual Arr_not_all_sides_oblivious_tag{};
struct Arr_has_open_side_tag :
    public virtual Arr_not_all_sides_oblivious_tag{};

struct Arr_all_sides_open_tag : public virtual Arr_not_all_sides_oblivious_tag{};

struct Arr_all_sides_not_open_tag {};
struct Arr_not_all_sides_not_open_tag {};

/* The following tags inherit from the Arr_not_all_sides_oblivious_tag;
 * They match any parameter space where at least one side of the boundary is not
 * oblivious.
 * The former matches any parameter space that none of its sides are either
 * identified, contracted, or closed. In other words, it matches parameter
 * spaces that cannot contain points on their boundaries.
 * The latter matches any parameter space that has at least one side that is
 * either an identified, contracted, or closed side. In other words, it matches
 * any parameter space, the boundary of which may contain a point.
 */
struct Arr_all_sides_not_finite_tag :
  public virtual Arr_not_all_sides_oblivious_tag {};
struct Arr_not_all_sides_not_finite_tag :
  public virtual Arr_not_all_sides_oblivious_tag {};

typedef std::true_type                                          Arr_true;
typedef std::false_type                                         Arr_false;

template <typename ArrSideCategory>
struct Arr_is_side_oblivious {
  typedef ArrSideCategory                                       Side_cat;
  typedef std::is_same<Side_cat, Arr_oblivious_side_tag>        Is_same;
  typedef std::bool_constant<Is_same::value>                    result;
  typedef typename result::type                                 type;
};

template <typename ArrSideCategory>
struct Arr_is_side_open {
  typedef ArrSideCategory                                       Side_cat;
  typedef std::is_same<Side_cat, Arr_open_side_tag>             Is_same;
  typedef std::bool_constant<Is_same::value>                    result;
  typedef typename result::type                                 type;
};

template <typename ArrSideCategory>
struct Arr_is_side_identified {
  typedef ArrSideCategory                                       Side_cat;
  typedef std::is_same<Side_cat, Arr_identified_side_tag>       Is_same;
  typedef std::bool_constant<Is_same::value>                    result;
  typedef typename result::type                                 type;
};

template <typename ArrSideCategory>
inline constexpr bool Arr_is_side_identified_v =
  Arr_is_side_identified<ArrSideCategory>::type::value;

template <typename ArrSideCategory>
struct Arr_is_side_contracted {
  typedef ArrSideCategory                                       Side_cat;
  typedef std::is_same<Side_cat, Arr_contracted_side_tag>       Is_same;
  typedef std::bool_constant<Is_same::value>                    result;
  typedef typename result::type                                 type;
};

template <typename ArrSideCategory>
struct Arr_is_side_closed {
  typedef ArrSideCategory                                       Side_cat;
  typedef std::is_same<Side_cat, Arr_closed_side_tag>           Is_same;
  typedef std::bool_constant<Is_same::value>                    result;
  typedef typename result::type                                 type;
};

/*! Struct to determine whether all side tags are "oblivious"
 */
template <typename ArrLeftSideCategory, typename ArrBottomSideCategory,
          typename ArrTopSideCategory, typename ArrRightSideCategory>
struct Arr_all_sides_oblivious_category {
  typedef ArrLeftSideCategory           Lef_side_cat;
  typedef ArrRightSideCategory          Rig_side_cat;
  typedef ArrBottomSideCategory         Bot_side_cat;
  typedef ArrTopSideCategory            Top_side_cat;

  typedef typename Arr_is_side_oblivious<Lef_side_cat>::result  Lef_obl;
  typedef typename Arr_is_side_oblivious<Rig_side_cat>::result  Rig_obl;
  typedef typename Arr_is_side_oblivious<Bot_side_cat>::result  Bot_obl;
  typedef typename Arr_is_side_oblivious<Top_side_cat>::result  Top_obl;

  /*! Boolean tag that is Arr_all_sides_oblivious_tag if all sides are
   * oblivious, otherwise Arr_not_all_sides_oblivious_tag
   */
  typedef std::conditional_t<Lef_obl::value && Rig_obl::value &&
                             Bot_obl::value && Top_obl::value,
                             Arr_all_sides_oblivious_tag,
                             Arr_not_all_sides_oblivious_tag>
    result;
};

/*! Struct to determine whether all side tags are "not-open"
 */
template <typename ArrLeftSideCategory, typename ArrBottomSideCategory,
          typename ArrTopSideCategory, typename ArrRightSideCategory>
struct Arr_all_sides_not_open_category {
public:
  typedef ArrLeftSideCategory           Lef_side_cat;
  typedef ArrRightSideCategory          Rig_side_cat;
  typedef ArrBottomSideCategory         Bot_side_cat;
  typedef ArrTopSideCategory            Top_side_cat;

private:
  typedef typename Arr_is_side_open<Lef_side_cat>::result       Lef_ope;
  typedef typename Arr_is_side_open<Rig_side_cat>::result       Rig_ope;
  typedef typename Arr_is_side_open<Bot_side_cat>::result       Bot_ope;
  typedef typename Arr_is_side_open<Top_side_cat>::result       Top_ope;

  static inline constexpr bool lef_not_ope = !Lef_ope::value;
  static inline constexpr bool rig_not_ope = !Rig_ope::value;
  static inline constexpr bool bot_not_ope = !Bot_ope::value;
  static inline constexpr bool top_not_ope = !Top_ope::value;

public:
  /*! Boolean tag that is Arr_all_sides_not_open_tag if all sides are not-open,
   * otherwise Arr_not_all_sides_not_open_tag
   */
  typedef std::conditional_t<lef_not_ope && rig_not_ope &&
                             bot_not_ope && top_not_ope,
                             Arr_all_sides_not_open_tag,
                             Arr_not_all_sides_not_open_tag>
    result;
};

/*! Struct to determine one of the following:
 * All sides are oblivious,
 * not all sides are oblivious, and not all sides are not open, or
 * Not all sides are oblivious, and all sides are not open
 * One distinct property of parameter spaces that are either oblivious or open
 * is that points cannot exist on the boundary of such parameter spaces.
 */
template <typename ArrLeftSideCategory, typename ArrBottomSideCategory,
          typename ArrTopSideCategory, typename ArrRightSideCategory>
struct Arr_sides_category {
public:
  typedef ArrLeftSideCategory                                   Lef_side_cat;
  typedef ArrBottomSideCategory                                 Bot_side_cat;
  typedef ArrTopSideCategory                                    Top_side_cat;
  typedef ArrRightSideCategory                                  Rig_side_cat;

private:
  typedef typename Arr_is_side_oblivious<Lef_side_cat>::result  Lef_obl;
  typedef typename Arr_is_side_oblivious<Rig_side_cat>::result  Rig_obl;
  typedef typename Arr_is_side_oblivious<Bot_side_cat>::result  Bot_obl;
  typedef typename Arr_is_side_oblivious<Top_side_cat>::result  Top_obl;

  typedef typename Arr_is_side_open<Lef_side_cat>::result       Lef_ope;
  typedef typename Arr_is_side_open<Rig_side_cat>::result       Rig_ope;
  typedef typename Arr_is_side_open<Bot_side_cat>::result       Bot_ope;
  typedef typename Arr_is_side_open<Top_side_cat>::result       Top_ope;

  static inline constexpr bool lef_obl_or_ope = Lef_obl::value || Lef_ope::value;
  static inline constexpr bool rig_obl_or_ope = Rig_obl::value || Rig_ope::value;
  static inline constexpr bool bot_obl_or_ope = Bot_obl::value || Bot_ope::value;
  static inline constexpr bool top_obl_or_ope = Top_obl::value || Top_ope::value;

  typedef std::conditional_t<lef_obl_or_ope &&
                             rig_obl_or_ope &&
                             bot_obl_or_ope &&
                             top_obl_or_ope,
                             Arr_all_sides_not_finite_tag,
                             Arr_not_all_sides_not_finite_tag>
    tmp;

public:
  typedef std::conditional_t<Lef_obl::value && Rig_obl::value &&
                             Bot_obl::value && Top_obl::value,
                             Arr_all_sides_oblivious_tag, tmp>
    result;
};

/*! Struct to check consistent tagging of identifications
 */
template <typename ArrLeftSideCategory, typename ArrBottomSideCategory,
          typename ArrTopSideCategory, typename ArrRightSideCategory>
struct Arr_sane_identified_tagging {
  typedef ArrLeftSideCategory                           Lef_side_cat;
  typedef ArrRightSideCategory                          Rig_side_cat;
  typedef ArrBottomSideCategory                         Bot_side_cat;
  typedef ArrTopSideCategory                            Top_side_cat;

  static inline constexpr bool lef_ide = Arr_is_side_identified_v<Lef_side_cat>;
  static inline constexpr bool rig_ide = Arr_is_side_identified_v<Rig_side_cat>;
  static inline constexpr bool bot_ide = Arr_is_side_identified_v<Bot_side_cat>;
  static inline constexpr bool top_ide = Arr_is_side_identified_v<Top_side_cat>;

  static inline constexpr bool lr_ide = lef_ide && rig_ide;
  static inline constexpr bool bt_ide = bot_ide && top_ide;

  static inline constexpr bool lr_not_ide = !lef_ide && !rig_ide;

  static inline constexpr bool bt_not_ide = !bot_ide && !top_ide;

  static inline constexpr bool lr_ok = lr_ide || lr_not_ide;
  static inline constexpr bool bt_ok = bt_ide || bt_not_ide;

  /*! Boolean tag that is bool_constant<true> if opposite sides are either
   * both identified or both not-identified,
   * otherwise bool_constant<false>
   */
  typedef std::bool_constant<lr_ok && bt_ok>              result;
  static inline constexpr bool value = result::value;
};

/*! Checks whether one of two boundary sides are identified
 * Observe that if one side is identified, the opposite side must be identified
 * as well. Thus:
 * (i)  When Arr_has_identified_sides is used to check whether two opposite
 *      sides are identified, the check for the second side is redundant.
 * (ii) When Arr_has_identified_sides is used to check whether two not-opposite
 *      sides are identified, the check applies to all four sides.
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_has_identified_sides {
  typedef ArrSideOneCategory            Side_one_cat;
  typedef ArrSideTwoCategory            Side_two_cat;

  typedef typename Arr_is_side_identified<Side_one_cat>::result Side_one_ide;
  typedef typename Arr_is_side_identified<Side_two_cat>::result Side_two_ide;

  /*! Boolean tag that is bool_constant<true> if one side is identified,
   * otherwise bool_constant<false>
   */
  typedef std::bool_constant<Side_one_ide::value || Side_two_ide::value> result;
};

/*! Checks whether one of two boundary sides are contracted
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_has_contracted_sides_two {
  typedef ArrSideOneCategory                                    Side_one_cat;
  typedef ArrSideTwoCategory                                    Side_two_cat;

  typedef typename Arr_is_side_contracted<Side_one_cat>::result Side_one_con;
  typedef typename Arr_is_side_contracted<Side_two_cat>::result Side_two_con;

  /*!\ Boolean tag that is bool_constant<true> if one side is identified,
   * otherwise bool_constant<false>
   */
  typedef std::bool_constant<Side_one_con::value ||
                             Side_two_con::value>               result;
};

/*! Checks whether one of two boundary sides are closed
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_has_closed_sides_two {
  typedef ArrSideOneCategory                                    Side_one_cat;
  typedef ArrSideTwoCategory                                    Side_two_cat;

  typedef typename Arr_is_side_closed<Side_one_cat>::result     Side_one_clo;
  typedef typename Arr_is_side_closed<Side_two_cat>::result     Side_two_clo;

  /*! Boolean tag that is bool_constant<true> if one side is identified,
   * otherwise bool_constant<false>
   */
  typedef std::bool_constant<Side_one_clo::value ||
                             Side_two_clo::value>               result;
};

/*! Checks whether one of two boundary sides are open
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_has_open_sides_two {
  typedef ArrSideOneCategory                                    Side_one_cat;
  typedef ArrSideTwoCategory                                    Side_two_cat;

  typedef typename Arr_is_side_open<Side_one_cat>::result       Side_one_ope;
  typedef typename Arr_is_side_open<Side_two_cat>::result       Side_two_ope;

  /*! Boolean tag that is bool_constant<true> if one side is identified,
   * otherwise bool_constant<false>
   */
  typedef std::bool_constant<Side_one_ope::value ||
                             Side_two_ope::value>               result;
};

/*! Categorizes two boundary sides:
 * If one side is identified           => Arr_has_identified_side_tag
 * Otherwise if one side is contracted => Arr_has_contracted_side_tag
 * Otherwise if one side is closed     => Arr_has_closed_side_tag
 * Otherwise if one side is open       => Arr_has_open_side_tag
 * Otherwise                           => Arr_all_sides_oblivious_tag
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_two_sides_category {
  typedef ArrSideOneCategory                            Side_one_cat;
  typedef ArrSideTwoCategory                            Side_two_cat;

  // One of the two sides is identified
  typedef typename Arr_has_identified_sides<Side_one_cat, Side_two_cat>::result
    Is_identified;

  // One of the two sides is contracted
  typedef typename Arr_has_contracted_sides_two<Side_one_cat,
                                                Side_two_cat>::result
    Is_contracted;

  // One of the two sides is closed
  typedef typename Arr_has_closed_sides_two<Side_one_cat, Side_two_cat>::result
    Is_closed;

  // One of the two sides is open
  typedef typename Arr_has_open_sides_two<Side_one_cat, Side_two_cat>::result
    Is_open;

public:
  typedef std::conditional_t<Is_identified::value, Arr_has_identified_side_tag,
    std::conditional_t<Is_contracted::value, Arr_has_contracted_side_tag,
      std::conditional_t<Is_closed::value, Arr_has_closed_side_tag,
        std::conditional_t<Is_open::value, Arr_has_open_side_tag,
          Arr_all_sides_oblivious_tag>>>>
    result;
};

/*! Categorizes all sides:
 * If one side is identified           => Arr_has_identified_side_tag
 * Otherwise if one side is contracted => Arr_has_contracted_side_tag
 * Otherwise if one side is closed     => Arr_has_closed_side_tag
 * Otherwise if one side is open       => Arr_has_open_side_tag
 * Otherwise (all sides oblivious)     => Arr_all_sides_oblivious_tag
 */
template <typename ArrSideOneCategory, typename ArrSideTwoCategory>
struct Arr_all_sides_category {};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
