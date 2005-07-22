// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_TAGS_H
#define CGAL_VORONOI_DIAGRAM_2_TAGS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/tags.h>

CGAL_BEGIN_NAMESPACE

//CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//====================================================================
//====================================================================

Tag_true operator&&(const Tag_true&, const Tag_true&) {
  return Tag_true();
}

Tag_false operator&&(const Tag_false&, const Tag_false&) {
  return Tag_false();
}

Tag_false operator&&(const Tag_true&, const Tag_false&) {
  return Tag_false();
}

Tag_false operator&&(const Tag_false&, const Tag_true&) {
  return Tag_false();
}


Tag_false operator||(const Tag_false&, const Tag_false&) {
  return Tag_false();
}

Tag_true operator||(const Tag_true&, const Tag_true&) {
  return Tag_true();
}

Tag_true operator||(const Tag_true&, const Tag_false&) {
  return Tag_true();
}

Tag_true operator||(const Tag_false&, const Tag_true&) {
  return Tag_true();
}

Tag_true  operator!(const Tag_false&) { return Tag_true(); }
Tag_false operator!(const Tag_true&) { return Tag_false(); }


//====================================================================
//====================================================================

template<bool b>
struct Boolean_tag {
  static const bool value = b;
};


template<bool b1, bool b2>
Boolean_tag<b1 && b2> operator&&(const Boolean_tag<b1>&,
				 const Boolean_tag<b2>&) {
  return Boolean_tag<b1 && b2>();
}


template<bool b1, bool b2>
Boolean_tag<b1 || b2> operator||(const Boolean_tag<b1>&,
				 const Boolean_tag<b2>&) {
  return Boolean_tag<b1 || b2>();
}

template<bool b>
Boolean_tag<!b> operator!(const Boolean_tag<b>&) { return Boolean_tag<!b>(); }


template<class Tag_t> struct To_boolean_tag;

template<>
struct To_boolean_tag<Tag_true> {
  typedef Boolean_tag<true> Tag;
};

template<>
struct To_boolean_tag<Tag_false> {
  typedef Boolean_tag<false> Tag;
};

template<class Tag_t>
struct Not {
  typedef Boolean_tag<!Tag_t::value> Tag;
};

template<class Tag1_t, class Tag2_t>
struct And {
  typedef Boolean_tag<Tag1_t::value && Tag2_t::value> Tag;
};

template<class Tag1_t, class Tag2_t>
struct Or {
  typedef Boolean_tag<Tag1_t::value || Tag2_t::value> Tag;
};

struct Tag_converter
{
  Boolean_tag<true> operator()(const Tag_true&) const {
    return Boolean_tag<true>();
  }

  Boolean_tag<false> operator()(const Tag_false&) const {
    return Boolean_tag<false>();
  }
};


//CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_TAGS_H
