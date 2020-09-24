// Copyright (c) 2020  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot

// This file is intentionally not protected against re-inclusion.
// It's aimed at being included from within a user code to
// make any structure inheriting from a face graph model a face graph
// model itself

// It is the responsibility of the including file to correctly set the
// macros CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME
// and optionally CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS.
// They are #undefed at the end of this file.

/**
* \ingroup PkgBGLRef
* \file CGAL/boost/graph/graph_traits_inheritance_macros.h
* Convenience header file defining the necessary specializations and overloads to make a
* class, inheriting from a model of a face graph concept, a model of that face graph concept itself.
* Prior to the inclusion of this header, specific macros must be defined and those macros will be
* undefined automatically when processing to the inclusion of this header.
* It is possible to include the header several times if the operation must be done for several classes.
* The macros that must be defined are the following:
*   - `CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME`: the inheriting class. If it is a template class, it must be instantiated parameters named as in `CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS` or parameters available in the scope including the header;
*   - `CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME`: the base class. it must be instantiated parameters named as in `CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS` or parameters available in the scope including the header;
*   - `CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS`: (optional) if the inheriting class, a list of template parameters separated by commas (`,`) including `class/typename/integral type`.
*
* Some examples are provided in \ref Surface_mesh/sm_derivation.cpp and \ref Polyhedron/poly_derivation.cpp.
*
*/

#include <CGAL/config.h>

#if !defined(CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME) || !defined(CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME)

CGAL_pragma_warning("\nBoth macros CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME and CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME "
                    "must be defined if you want to use this file\n")

#else

#ifdef CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS

namespace boost {

template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS>
struct graph_traits<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME> :
  public graph_traits<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME>
{};

template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS, typename CGAL_XX_YATP>
struct property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, CGAL_XX_YATP> :
  public property_map<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, CGAL_XX_YATP>
{};

} // boost namespace

#define CGAL_PM_DT_SPEC(DTAG) \
namespace boost {\
template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS, typename CGAL_XX_YATP> \
struct property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> > \
  : property_map<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, DTAG<CGAL_XX_YATP> > \
{};\
} /* boost namespace */\
\
namespace CGAL { \
template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS, typename CGAL_XX_YATP>\
typename boost::property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> >::type \
get(DTAG<CGAL_XX_YATP> t, CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME& g) \
{ \
  return get(t, static_cast<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME&>(g)); \
} \
\
template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS, typename CGAL_XX_YATP>\
typename boost::property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> >::const_type \
get(DTAG<CGAL_XX_YATP> t, const CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME& g) \
{ \
  return get(t, static_cast<const CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME&>(g)); \
}\
} //CGAL namespace

CGAL_PM_DT_SPEC(CGAL::dynamic_vertex_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_halfedge_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_face_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_edge_property_t)

#undef CGAL_PM_DT_SPEC

namespace CGAL {

template <CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS, typename CGAL_XX_YATP>
struct graph_has_property<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, CGAL_XX_YATP> :
  public CGAL::graph_has_property<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, CGAL_XX_YATP>
{};


} // CGAL namespace

#undef CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS

#else

namespace boost {

template <>
struct graph_traits<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME> :
  public graph_traits<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME>
{};

template <typename CGAL_XX_YATP>
struct property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, CGAL_XX_YATP> :
  public property_map<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, CGAL_XX_YATP>
{};

} // boost namespace

#define CGAL_PM_DT_SPEC(DTAG) \
namespace boost {\
template <typename CGAL_XX_YATP> \
struct property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> > \
  : property_map<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, DTAG<CGAL_XX_YATP> > \
{};\
} /* boost namespace */\
\
namespace CGAL { \
template <typename CGAL_XX_YATP>\
typename boost::property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> >::type \
get(DTAG<CGAL_XX_YATP> t, CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME& g) \
{ \
  return get(t, static_cast<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME&>(g)); \
} \
\
template <typename CGAL_XX_YATP>\
typename boost::property_map<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, DTAG<CGAL_XX_YATP> >::const_type \
get(DTAG<CGAL_XX_YATP> t, const CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME& g) \
{ \
  return get(t, static_cast<const CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME&>(g)); \
}\
} //CGAL namespace

CGAL_PM_DT_SPEC(CGAL::dynamic_vertex_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_halfedge_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_face_property_t)
CGAL_PM_DT_SPEC(CGAL::dynamic_edge_property_t)

#undef CGAL_PM_DT_SPEC

namespace CGAL {

template <typename CGAL_XX_YATP>
struct graph_has_property<CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME, CGAL_XX_YATP> :
  public CGAL::graph_has_property<CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME, CGAL_XX_YATP>
{};


} // CGAL namespace

#endif


#undef CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME
#undef CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME

#endif

