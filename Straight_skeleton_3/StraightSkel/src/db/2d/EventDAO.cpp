/**
 * @file   db/2d/EventDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-04
 */

#include "db/2d/EventDAO.h"

#include "data/2d/skel/ConstOffsetEvent.h"
#include "data/2d/skel/EdgeEvent.h"
#include "data/2d/skel/SplitEvent.h"
#include "data/2d/skel/TriangleEvent.h"
#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/2d/NodeDAO.h"
#include "db/2d/StraightSkeletonDAO.h"

namespace db { namespace _2d {

EventDAO::EventDAO() {
    // intentionally does nothing
}

EventDAO::~EventDAO() {
    // intentionally does nothing
}

std::string EventDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Events (\n"
            "  SkelID INTEGER NOT NULL,\n"
            "  EventID INTEGER NOT NULL,\n"
            "  etype INTEGER,\n"
            "  NID INTEGER,\n"
            "  PRIMARY KEY (SkelID, EventID)\n"
            ");");
    return schema;
}

int EventDAO::nextEventID(int skelid) {
    int eventid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(EventID) FROM Events WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        eventid = 1;
        if (stmt->execute()) {
            eventid = stmt->getInteger(0) + 1;
        }
    }
    return eventid;
}

int EventDAO::insert(AbstractEventSPtr event) {
    int result = -1;
    if (!event->getSkel()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int skelid = event->getSkel()->getID();
    if (skelid <= 0) {
        StraightSkeletonDAOSPtr dao_skel = DAOFactory::getStraightSkeletonDAO();
        skelid = dao_skel->createSkelID(event->getSkel());
    }
    if (skelid > 0) {
        NodeSPtr node;
        if (event->getType() == AbstractEvent::EDGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_2d::skel::EdgeEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::SPLIT_EVENT) {
            node = std::dynamic_pointer_cast<data::_2d::skel::SplitEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
            node = std::dynamic_pointer_cast<data::_2d::skel::TriangleEvent>(event)->getNode();
        }
        int eventid = nextEventID(skelid);
        if (node) {
            NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
            if (node->getID() > 0) {
                dao_node->update(node);
            } else {
                dao_node->insert(node);
            }
            std::string sql("INSERT INTO Events (SkelID, EventID, etype, NID) "
                    "VALUES (?, ?, ?, ?);");
            SQLiteStmtSPtr stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->bindInteger(2, eventid);
                stmt->bindInteger(3, event->getType());
                stmt->bindInteger(4, node->getID());
                if (stmt->execute() > 0) {
                    event->setID(eventid);
                    result = eventid;
                }
            }
        } else {
            std::string sql("INSERT INTO Events (SkelID, EventID, etype) "
                    "VALUES (?, ?, ?);");
            SQLiteStmtSPtr stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->bindInteger(2, eventid);
                stmt->bindInteger(3, event->getType());
                if (stmt->execute() > 0) {
                    event->setID(eventid);
                    result = eventid;
                }
            }
        }
    }
    return result;
}

bool EventDAO::del(AbstractEventSPtr event) {
    bool result = false;
    if (!event->getSkel()) {
        return false;
    }
    int skelid = event->getSkel()->getID();
    int eventid = event->getID();
    if (skelid < 0 || eventid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Events WHERE SkelID=? AND EventID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, eventid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    if (result) {
        event->setID(-1);
    }
    return result;
}

AbstractEventSPtr EventDAO::find(int skelid, int eventid) {
    AbstractEventSPtr result;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT etype, NID FROM Events WHERE SkelID=? AND EventID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, eventid);
        if (stmt->execute() > 0) {
            int etype = stmt->getInteger(0);
            int nid = stmt->getInteger(1);
            NodeSPtr node;
            if (nid > 0) {
                NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
                node = dao_node->find(skelid, nid);
            }
            if (etype == AbstractEvent::CONST_OFFSET_EVENT) {
                result = data::_2d::skel::ConstOffsetEvent::create();
            } else if (etype == AbstractEvent::EDGE_EVENT) {
                data::_2d::skel::EdgeEventSPtr edge_event =
                        data::_2d::skel::EdgeEvent::create();
                edge_event->setNode(node);
                result = edge_event;
            } else if (etype == AbstractEvent::SPLIT_EVENT) {
                data::_2d::skel::SplitEventSPtr split_event =
                        data::_2d::skel::SplitEvent::create();
                split_event->setNode(node);
                result = split_event;
            } else if (etype == AbstractEvent::TRIANGLE_EVENT) {
                data::_2d::skel::TriangleEventSPtr triangle_event =
                        data::_2d::skel::TriangleEvent::create();
                triangle_event->setNode(node);
                result = triangle_event;
            }
            result->setID(eventid);
        }
    }
    return result;
}

bool EventDAO::update(AbstractEventSPtr event) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
