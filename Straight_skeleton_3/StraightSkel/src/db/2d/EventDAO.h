/**
 * @file   db/2d/EventDAO.h
 * @author Gernot Walzl
 * @date   2013-06-04
 */

#ifndef DB_2D_EVENTDAO_H
#define DB_2D_EVENTDAO_H

#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/AbstractEvent.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::skel::AbstractEvent;
using data::_2d::skel::AbstractEventSPtr;

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

#endif /* DB_2D_EVENTDAO_H */
