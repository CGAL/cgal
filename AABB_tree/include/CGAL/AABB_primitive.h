// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_AABB_PRIMITIVE_H
#define CGAL_AABB_PRIMITIVE_H

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/AABB_tree/internal/Has_nested_type_Shared_data.h>
#include <CGAL/property_map.h>
#include <CGAL/tags.h>

namespace CGAL {

//class for the typedefs
template < class Id_,
           class ObjectPropertyMap,
           class PointPropertyMap >
struct AABB_primitive_base
{
  typedef typename boost::property_traits< ObjectPropertyMap >::value_type Datum; //datum type
  typedef typename boost::property_traits< PointPropertyMap  >::value_type Point; //point type
  typedef typename boost::property_traits< ObjectPropertyMap >::reference Datum_reference; //reference datum type
  typedef typename boost::property_traits< PointPropertyMap  >::reference Point_reference; //reference point type
  typedef Id_ Id; // Id type

protected:
  Id m_id;

public:
  // constructors
  AABB_primitive_base(Id id) : m_id(id) {}

  Id id() const {return m_id;}
};

#ifdef DOXYGEN_RUNNING
/*!
 * \ingroup PkgAABBTreeRef
 * Generic primitive type.
 * The two property maps which are template parameters of the class enable to get the datum and the reference point of
 * the primitive from the identifier. The last template parameter controls whether the primitive class holds a copy of the datum.
 *
 * \cgalModels `AABBPrimitive` if `ExternalPropertyMaps` is `CGAL::Tag_false`.
 * \cgalModels `AABBPrimitiveWithSharedData` if `ExternalPropertyMaps` is `CGAL::Tag_true`.
 *
 * \tparam ObjectPropertyMap is a model of `ReadablePropertyMap` with `Id` as
 *           `key_type`. It must be a model of `CopyConstructible`, `DefaultConstructible`, and `CopyAssignable`.
 * \tparam PointPropertyMap is a model of `ReadablePropertyMap` with `Id` as
 *           `key_type`. It must be a model of `CopyConstructible`, `DefaultConstructible`, and `CopyAssignable`.
 * \tparam ExternalPropertyMaps either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former
 *          case, the property maps will be stored in the traits class, while
 *          in the latter they will be stored in the primitive
 *          (which increases the size of each primitive).
 * \tparam CacheDatum either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case,
 *           the datum is stored in the primitive, while in the latter
 *           it is constructed on the fly to reduce the memory footprint.
 *           The default is `CGAL::Tag_false` (datum is not stored).
 *
 * \sa `AABB_segment_primitive<Iterator,CacheDatum>`
 * \sa `AABB_triangle_primitive<Iterator,CacheDatum>`
 * \sa `AABB_halfedge_graph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatum>`
 * \sa `AABB_face_graph_triangle_primitive<FaceGraph,OneFaceGraphPerTree,CacheDatum>`
 */
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap,
            class ExternalPropertyMaps,
            class CacheDatum>
struct AABB_primitive
{
  /// \name Types
  /// @{
  /*!
  The datum type.
  */
  typedef boost::property_traits< ObjectPropertyMap >::value_type Datum;

  /*!
  The point type.
  */
  typedef boost::property_traits< PointPropertyMap >::value_type Point;

  /*!
  The reference datum type.
  */
  typedef boost::property_traits< ObjectPropertyMap >::reference Datum_reference;

  /*!
  The reference point type.
  */
  typedef boost::property_traits< PointPropertyMap >::reference Point_reference;

  /*!
  Id type.
  */
  typedef Id Id;
  /// @}

  /*!
  constructs a primitive and initializes the property maps.
  */
  AABB_primitive(Id id,
                 ObjectPropertyMap o_pmap=ObjectPropertyMap(),
                 PointPropertyMap p_pmap=PointPropertyMap());

  /*!
  constructs a primitive from an iterator with `Id` as value type
  and initializes the property maps.
  */
  template <class Iterator>
  AABB_primitive(Iterator it,
                 ObjectPropertyMap o_pmap=ObjectPropertyMap(),
                 PointPropertyMap p_pmap=PointPropertyMap());
};
#else
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap,
            class ExternalPropertyMaps,
            class CacheDatum>
class AABB_primitive;


//no caching, property maps internally stored
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap >
class AABB_primitive<Id, ObjectPropertyMap, PointPropertyMap,Tag_false,Tag_false>
  : public AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap>
{
  typedef AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap> Base;
  ObjectPropertyMap m_obj_pmap;
  PointPropertyMap m_pt_pmap;
public:
  AABB_primitive(Id id, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap pt_pmap=PointPropertyMap())
    : Base(id), m_obj_pmap(obj_pmap), m_pt_pmap(pt_pmap) {}

  template <class Iterator>
  AABB_primitive(Iterator it, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap pt_pmap=PointPropertyMap())
    : Base(*it), m_obj_pmap(obj_pmap), m_pt_pmap(pt_pmap) {}

  typename Base::Datum_reference
  datum() const { return get(m_obj_pmap,this->m_id); }

  typename Base::Point_reference
  reference_point() const { return get(m_pt_pmap,this->m_id); }
};

//caching, property maps internally stored
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap >
class AABB_primitive<Id, ObjectPropertyMap, PointPropertyMap,Tag_false,Tag_true>
  : public AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap>
{
  typedef AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap> Base;
  typename boost::property_traits< ObjectPropertyMap >::value_type m_datum;
  PointPropertyMap m_pt_pmap;
public:
  typedef const typename Base::Datum& Datum_reference;

  AABB_primitive(Id id, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap pt_pmap=PointPropertyMap())
    : Base(id), m_datum( get(obj_pmap,id) ), m_pt_pmap(pt_pmap){}

  template <class Iterator>
  AABB_primitive(Iterator it, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap pt_pmap=PointPropertyMap())
    : Base(*it), m_datum( get(obj_pmap,*it) ), m_pt_pmap(pt_pmap){}


  Datum_reference datum() const { return m_datum; }

  typename Base::Point_reference
  reference_point() const { return get(m_pt_pmap,this->m_id); }
};

//no caching, property maps are stored outside the class
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap >
class AABB_primitive<Id, ObjectPropertyMap, PointPropertyMap,Tag_true,Tag_false>
  : public AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap>
{
  typedef AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap> Base;
public:
  typedef std::pair<ObjectPropertyMap,PointPropertyMap> Shared_data;

  AABB_primitive(Id id, ObjectPropertyMap=ObjectPropertyMap(), PointPropertyMap=PointPropertyMap())
    : Base(id) {}

  template <class Iterator>
  AABB_primitive(Iterator it, ObjectPropertyMap=ObjectPropertyMap(), PointPropertyMap=PointPropertyMap())
    : Base(*it) {}

  typename Base::Datum_reference
  datum(const Shared_data& data) const { return get(data.first,this->m_id); }

  typename Base::Point_reference
  reference_point(const Shared_data& data) const { return get(data.second,this->m_id); }

  static Shared_data construct_shared_data(ObjectPropertyMap obj, PointPropertyMap pt) {return Shared_data(obj,pt);}
};


//caching, property map is stored outside the class
template <  class Id,
            class ObjectPropertyMap,
            class PointPropertyMap >
class AABB_primitive<Id, ObjectPropertyMap, PointPropertyMap,Tag_true,Tag_true>
  : public AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap>
{
  typedef AABB_primitive_base<Id,ObjectPropertyMap,PointPropertyMap> Base;
  typename boost::property_traits< ObjectPropertyMap >::value_type m_datum;
public:
  typedef PointPropertyMap Shared_data;
  typedef const typename Base::Datum& Datum_reference;

  AABB_primitive(Id id, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap=PointPropertyMap())
    : Base(id), m_datum( get(obj_pmap,id) ) {}

  template <class Iterator>
  AABB_primitive(Iterator it, ObjectPropertyMap obj_pmap=ObjectPropertyMap(), PointPropertyMap=PointPropertyMap())
    : Base(*it), m_datum( get(obj_pmap,*it) ) {}

  Datum_reference datum(Shared_data) const { return m_datum; }

  typename Base::Point_reference
  reference_point(const Shared_data& data) const { return get(data,this->m_id); }

  static Shared_data construct_shared_data(ObjectPropertyMap, PointPropertyMap pt) {return pt;}
};

#endif //DOXYGEN_RUNNING

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_PRIMITIVE_H

