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
 * @file   data/3d/skel/SphericalInversionEvent.h
 * @author Gernot Walzl
 * @date   2013-03-13
 */

#ifndef DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H
#define DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalInversionEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalInversionEvent();
    static SphericalInversionEventSPtr create();
    CGAL::FT getOffset() const;
    void setOffset(CGAL::FT offset);
    SphericalPolygonSPtr getPolygon() const;
    void setPolygon(SphericalPolygonSPtr polygon);
protected:
    SphericalInversionEvent();
    CGAL::FT offset_;
    SphericalPolygonSPtr polygon_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALINVERSIONEVENT_H */

