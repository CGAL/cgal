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
 * @file   data/3d/skel/SkelFacetData.h
 * @author Gernot Walzl
 * @date   2102-05-07
 */

#ifndef DATA_3D_SKEL_SKELFACETDATA_H
#define DATA_3D_SKEL_SKELFACETDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/FacetData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SkelFacetData : public FacetData {
public:
    virtual ~SkelFacetData();

    static SkelFacetDataSPtr create(FacetSPtr facet);

    FacetSPtr getOffsetFacet() const;
    void setOffsetFacet(FacetSPtr offset_facet);

    FacetSPtr getFacetOrigin() const;
    void setFacetOrigin(FacetSPtr facet_origin);

    CGAL::FT getSpeed() const;
    void setSpeed(CGAL::FT speed);

    int getStepID() const;
    void setStepID(int id);

protected:
    SkelFacetData();
    FacetWPtr offset_facet_;
    FacetWPtr facet_origin_;
    CGAL::FT speed_;

    // This is the ID of the last event that was processed (popped from the queue)
    // which impacted this facet.
    // The point is that we pop the queue, we can compare the step ID of the facet involved
    // in the handling of the event to the step contained in the event.
    // If the step of the event is *smaller* than the step of the facet, it means the facet
    // has been modified after the event was created, so we ignore the event.
    int step_id_ = -1;
};

} } }

#endif /* DATA_3D_SKEL_SKELFACETDATA_H */
