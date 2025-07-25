// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   data/3d/skel/AbstractEvent.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_ABSTRACTEVENT_H
#define DATA_3D_SKEL_ABSTRACTEVENT_H

#include "debug.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/Vertex.h"
#include "data/3d/Edge.h"
#include "data/3d/Facet.h"

#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class AbstractEvent {
public:
    virtual ~AbstractEvent();

    PolyhedronSPtr getPolyhedronResult() const;
    void setPolyhedronResult(PolyhedronSPtr polyhedron);
    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<AbstractEventSPtr>::iterator getListIt() const;
    void setListIt(std::list<AbstractEventSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    int getStepID() const;
    void setStepID(int id);

    virtual void setHighlight(bool highlight);
    virtual CGAL::FT getOffset() const = 0; // abstract

    static const int SAVE_OFFSET_EVENT = 0;

    static const int CONST_OFFSET_EVENT = 1;

    /** generic vanish event */
    static const int VANISH_EVENT = -1;

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

    virtual int getType() const;

    virtual std::string toString() const;

    virtual bool isValid() const;

    virtual bool isObsolete() const;

protected:
    AbstractEvent();

    PolyhedronSPtr polyhedron_result_;
    StraightSkeletonWPtr skel_;
    std::list<AbstractEventSPtr>::iterator list_it_;
    int type_;
    int id_; // id of the event
    int step_id_ = -1; // id of the *step* at which the event was created

private:
    static int next_id_;
};

class AbstractEventSPtrCompare
  : public CGAL::cpp98::binary_function<bool, AbstractEventSPtr, AbstractEventSPtr>
{
public:
    bool operator()(const AbstractEventSPtr& eventA,
                    const AbstractEventSPtr& eventB) const
    {
      if (eventA->getOffset() == eventB->getOffset()) {
          if (eventA->getType() == eventB->getType()) { // @fixme get rid of that and IDs directly?
              // Give priority to newer (higher) IDs. The point is that if an event has been updated
              // to a different type and appears multiple (non-zombie) times, it will be processed
              // with the updated type.
              return eventA->getID() < eventB->getID();
          }
          return (eventA->getType() > eventB->getType());
      }

      return (eventA->getOffset() < eventB->getOffset());
    }
};

struct VertexFacetNeighborhood
{
    VertexFacetNeighborhood() {}
    VertexFacetNeighborhood(const VertexSPtr v)
        : incident_facets_(VertexFacetNeighborhood::collectIncidentFacets(v))
    { }

    static std::array<int, 3> collectIncidentFacets(const VertexSPtr vertex)
    {
        CGAL_SS3_DEBUG_SPTR(vertex);
        std::array<int, 3> result;
        const std::list<FacetWPtr>& ifs = vertex->facets();
        CGAL_assertion(ifs.size() == 3);
        std::list<FacetWPtr>::const_iterator it = ifs.begin();
        for (int i = 0; i < 3; ++i, ++it) {
            if (FacetSPtr f = it->lock()) {
                result[i] = f->getBasePlaneID();
            } else {
                CGAL_assertion(false);
                result[i] = -1;
            }
        }
        return result;
    }

    bool checkNeighborhoodConsistency(const VertexSPtr other) const
    {
        CGAL_SS3_DEBUG_SPTR(other);
        const std::array<int, 3>& other_ifs = collectIncidentFacets(other);
        const bool result = (incident_facets_ == other_ifs);
        // std::cout << "  Compare: " << incident_facets_[0] << " (old) " << other_ifs[0] << " (new)" << std::endl;
        // std::cout << "  Compare: " << incident_facets_[1] << " (old) " << other_ifs[1] << " (new)" << std::endl;
        // std::cout << "  Compare: " << incident_facets_[2] << " (old) " << other_ifs[2] << " (new)" << std::endl;
        if (!result) {
        }
        return result;
    }

private:
    std::array<int, 3> incident_facets_;
};

struct EdgeFacetNeighborhood
{
    EdgeFacetNeighborhood() {}
    EdgeFacetNeighborhood(const EdgeSPtr edge)
        : incident_facets_(collectIncidentFacets(edge))
    { }

    static std::array<int, 4> collectIncidentFacets(const EdgeSPtr& edge)
    {
        CGAL_SS3_DEBUG_SPTR(edge);
        std::array<int, 4> result = {{ edge->getFacetL()->getBasePlaneID(),
                                       edge->getFacetR()->getBasePlaneID(),
                                       edge->getFacetSrc()->getBasePlaneID(),
                                       edge->getFacetDst()->getBasePlaneID() }};
        return result;
    }

    bool checkNeighborhoodConsistency(const EdgeSPtr other) const
    {
        CGAL_SS3_DEBUG_SPTR(other);

        const std::array<int, 4>& other_ifs = collectIncidentFacets(other);
        const bool result = (incident_facets_ == other_ifs);
        if (!result) {
            // std::cout << "  Compare: " << incident_facets_[0] << " (old) " << other_ifs[0] << " (new)" << std::endl;
            // std::cout << "  Compare: " << incident_facets_[1] << " (old) " << other_ifs[1] << " (new)" << std::endl;
            // std::cout << "  Compare: " << incident_facets_[2] << " (old) " << other_ifs[2] << " (new)" << std::endl;
            // std::cout << "  Compare: " << incident_facets_[3] << " (old) " << other_ifs[3] << " (new)" << std::endl;
        }
        return result;
    }

private:
    std::array<int, 4> incident_facets_;
};

} } }

#endif /* DATA_3D_ABSTRACTEVENT_H */
