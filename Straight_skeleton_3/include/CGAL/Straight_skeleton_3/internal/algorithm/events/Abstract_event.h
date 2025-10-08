// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * file   data/3d/skel/AbstractEvent.h
 * author Gernot Walzl
 * date   2012-03-27
 */

#ifndef CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_EVENT_H
#define CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_EVENT_H

#include <CGAL/Straight_skeleton_3/internal/debug.h>
#include <CGAL/Straight_skeleton_3/IO/String_factory.h>

#include <CGAL/number_utils.h>

#include <iostream>
#include <memory>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace internal {
namespace algorithm {

template <typename Traits>
class AbstractEvent
{
  using AbstractEventSPtr = std::shared_ptr<AbstractEvent<Traits> >;

  using FT = typename Traits::FT;

public:
  AbstractEvent(int type = -2)
    : type_(type),
      id_(next_id_++)
  { }

  virtual ~AbstractEvent() { /*intentionally does nothing*/ }

  typename std::list<AbstractEventSPtr>::iterator getListIt() const
  {
    return this->list_it_;
  }

  void setListIt(typename std::list<AbstractEventSPtr>::iterator list_it)
  {
    this->list_it_ = list_it;
  }

  int getID() const
  {
    return this->id_;
  }

  void setID(const int id)
  {
    this->id_ = id;
  }

  const FT& getTime() const
  {
    return time_;
  }

  void setTime(const FT& time)
  {
    this->time_ = time;
  }

  static const int SAVE_EVENT = 0;

  static const int CONST_TIME_EVENT = 1;

  /** generic vanish event */
  static const int VANISH_EVENT = -1; // @tmp give it a proper ID

  /** 1 edge vanish event */
  static const int EDGE_EVENT = 2;

  /** 2 edge vanish event */
  static const int EDGE_MERGE_EVENT = 3;

  /** 3 edge vanish event */
  static const int TRIANGLE_EVENT = 4;

  /** 4 edge vanish event */
  static const int DBL_EDGE_MERGE_EVENT = 5;

  /** 5 edge vanish event */
  static const int DBL_TRIANGLE_EVENT = 6;

  /** 6 edge vanish event */
  static const int TETRAHEDRON_EVENT = 7;

  /** vertex-vertex contact event I */
  static const int VERTEX_EVENT = 8;

  /** vertex-vertex contact event II */
  static const int FLIP_VERTEX_EVENT = 9;

  /** vertex-edge contact event */
  static const int SURFACE_EVENT = 10;

  /** vertex-vertex-edge contact event I */
  static const int POLYHEDRON_SPLIT_EVENT = 11;

  /** vertex-vertex-edge contact event II */
  static const int SPLIT_MERGE_EVENT = 12;

  /** edge-edge contact event */
  static const int EDGE_SPLIT_EVENT = 13;

  /** vertex-facet contact event */
  static const int PIERCE_EVENT = 14;

  int getType() const {
    return this->type_;
  }

  virtual bool isValid() const
  {
    return true;
  }

  virtual bool isObsolete() const
  {
    return false;
  }

  virtual std::string toString() const
  {
    std::stringstream sstr;
    switch (getType()) {
      case CONST_TIME_EVENT:
        sstr << "ConstTimeEvent";
        break;
      case SAVE_EVENT:
        sstr << "SaveEvent";
        break;
      case VANISH_EVENT:
        sstr << "VanishEvent";
        break;
      case EDGE_EVENT:
        sstr << "EdgeEvent";
        break;
      case EDGE_MERGE_EVENT:
        sstr << "EdgeMergeEvent";
        break;
      case TRIANGLE_EVENT:
        sstr << "TriangleEvent";
        break;
      case DBL_EDGE_MERGE_EVENT:
        sstr << "DblEdgeMergeEvent";
        break;
      case DBL_TRIANGLE_EVENT:
        sstr << "DblTriangleEvent";
        break;
      case TETRAHEDRON_EVENT:
        sstr << "TetrahedronEvent";
        break;
      case VERTEX_EVENT:
        sstr << "VertexEvent";
        break;
      case FLIP_VERTEX_EVENT:
        sstr << "FlipVertexEvent";
        break;
      case SURFACE_EVENT:
        sstr << "SurfaceEvent";
        break;
      case POLYHEDRON_SPLIT_EVENT:
        sstr << "PolyhedronSplitEvent";
        break;
      case SPLIT_MERGE_EVENT:
        sstr << "SplitMergeEvent";
        break;
      case EDGE_SPLIT_EVENT:
        sstr << "EdgeSplitEvent";
        break;
      case PIERCE_EVENT:
        sstr << "PierceEvent";
        break;
      default:
        sstr << "AbstractEvent";
    }
    sstr << "(offset=" << IO::StringFactory::fromDouble(CGAL::to_double(getTime())) << ")";
    return sstr.str();
  }

protected:
  typename std::list<AbstractEventSPtr>::iterator list_it_;

  int type_;
  int id_; // id of the event
  FT time_;

private:
  static int next_id_;
};

template <typename Traits>
int AbstractEvent<Traits>::next_id_ = 0;

template <typename Traits>
class AbstractEventSPtrCompare
  : public CGAL::cpp98::binary_function<bool,
                                        std::shared_ptr<AbstractEvent<Traits> >,
                                        std::shared_ptr<AbstractEvent<Traits> > >
{
public:
  bool operator()(const std::shared_ptr<AbstractEvent<Traits> >& eventA,
                  const std::shared_ptr<AbstractEvent<Traits> >& eventB) const
  {
    CGAL_SS3_DEBUG_SPTR(eventA);
    CGAL_SS3_DEBUG_SPTR(eventB);
    if (eventA->getTime() == eventB->getTime()) {
      if (eventA->getType() == eventB->getType()) { // @fixme get rid of that and IDs directly?
        // Give priority to newer (higher) IDs. The point is that if an event has been updated
        // to a different type and appears multiple (non-zombie) times, it will be processed
        // with the updated type.
        return eventA->getID() < eventB->getID();
      }
      return (eventA->getType() > eventB->getType());
    }

    return (eventA->getTime() < eventB->getTime());
  }
};

template <typename Traits>
struct VertexFacetNeighborhood
{
  using Polyhedron = HDS::Polyhedron<Traits>;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using FacetWPtr = typename Polyhedron::FacetWPtr;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  VertexFacetNeighborhood() {}
  VertexFacetNeighborhood(const VertexSPtr& v)
    : incident_facets_(VertexFacetNeighborhood::collectIncidentFacets(v))
  { }

  static std::array<int, 3> collectIncidentFacets(const VertexSPtr& vertex)
  {
    CGAL_SS3_DEBUG_SPTR(vertex);
    std::array<int, 3> result;
    const std::list<FacetWPtr>& ifs = vertex->facets();
    CGAL_assertion(ifs.size() == 3);
    typename std::list<FacetWPtr>::const_iterator it = ifs.begin();
    for (int i = 0; i < 3; ++i, ++it) {
      if (FacetSPtr f = it->lock()) {
        result[i] = f->getID();
      } else {
        CGAL_assertion(false);
        result[i] = -1;
      }
    }
    return result;
  }

  bool checkNeighborhoodConsistency(const VertexSPtr& other) const
  {
    CGAL_SS3_DEBUG_SPTR(other);
    const std::array<int, 3>& other_ifs = collectIncidentFacets(other);
    const bool result = (incident_facets_ == other_ifs);
    // std::cout << "  Compare: " << incident_facets_[0] << " (old) " << other_ifs[0] << " (new)" << std::endl;
    // std::cout << "  Compare: " << incident_facets_[1] << " (old) " << other_ifs[1] << " (new)" << std::endl;
    // std::cout << "  Compare: " << incident_facets_[2] << " (old) " << other_ifs[2] << " (new)" << std::endl;
    return result;
  }

private:
  std::array<int, 3> incident_facets_;
};

template <typename Traits>
struct EdgeFacetNeighborhood
{
  using Polyhedron = HDS::Polyhedron<Traits>;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;

  EdgeFacetNeighborhood() {}
  EdgeFacetNeighborhood(const EdgeSPtr edge)
    : incident_facets_(collectIncidentFacets(edge))
  { }

  static std::array<int, 4> collectIncidentFacets(const EdgeSPtr& edge)
  {
    CGAL_SS3_DEBUG_SPTR(edge);
    std::array<int, 4> result = {{ edge->getFacetL()->getID(),
                                   edge->getFacetR()->getID(),
                                   edge->getFacetSrc()->getID(),
                                   edge->getFacetDst()->getID() }};
    return result;
  }

  bool checkNeighborhoodConsistency(const EdgeSPtr& other) const
  {
    CGAL_SS3_DEBUG_SPTR(other);

    const std::array<int, 4>& other_ifs = collectIncidentFacets(other);
    const bool result = (incident_facets_ == other_ifs);
    // if (!result) {
    //   std::cout << "  Compare: " << incident_facets_[0] << " (old) " << other_ifs[0] << " (new)" << std::endl;
    //   std::cout << "  Compare: " << incident_facets_[1] << " (old) " << other_ifs[1] << " (new)" << std::endl;
    //   std::cout << "  Compare: " << incident_facets_[2] << " (old) " << other_ifs[2] << " (new)" << std::endl;
    //   std::cout << "  Compare: " << incident_facets_[3] << " (old) " << other_ifs[3] << " (new)" << std::endl;
    // }
    return result;
  }

private:
  std::array<int, 4> incident_facets_;
};

} // namespace algorithm
} // namespace internal
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_INTERNAL_ALGORITHM_ABSTRACT_EVENT_H */
