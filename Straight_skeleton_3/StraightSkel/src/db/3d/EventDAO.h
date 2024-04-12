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
 * @file   db/3d/EventDAO.h
 * @author Gernot Walzl
 * @date   2013-06-04
 */

#ifndef DB_3D_EVENTDAO_H
#define DB_3D_EVENTDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::skel::AbstractEvent;
using data::_3d::skel::AbstractEventSPtr;

class EventDAO {
friend class DAOFactory;
public:
    virtual ~EventDAO();
    std::string getTableSchema() const;
    int insert(AbstractEventSPtr event);
    bool del(AbstractEventSPtr event);
    AbstractEventSPtr find(int skelid, int eventid);
    bool update(AbstractEventSPtr event);
protected:
    EventDAO();
    int nextEventID(int skelid);
};

} }

#endif /* DB_3D_EVENTDAO_H */
