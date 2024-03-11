/**
 * @file   db/3d/EventDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-04
 */

#include "db/3d/EventDAO.h"

#include "data/3d/skel/ConstOffsetEvent.h"
#include "data/3d/skel/EdgeEvent.h"
#include "data/3d/skel/EdgeMergeEvent.h"
#include "data/3d/skel/TriangleEvent.h"
#include "data/3d/skel/DblEdgeMergeEvent.h"
#include "data/3d/skel/DblTriangleEvent.h"
#include "data/3d/skel/TetrahedronEvent.h"
#include "data/3d/skel/VertexEvent.h"
#include "data/3d/skel/FlipVertexEvent.h"
#include "data/3d/skel/SurfaceEvent.h"
#include "data/3d/skel/PolyhedronSplitEvent.h"
#include "data/3d/skel/SplitMergeEvent.h"
#include "data/3d/skel/EdgeSplitEvent.h"
#include "data/3d/skel/PierceEvent.h"
#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/NodeDAO.h"
#include "db/3d/StraightSkeletonDAO.h"

namespace db { namespace _3d {

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
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeMergeEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::TriangleEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::DblEdgeMergeEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::DblTriangleEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::TetrahedronEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::VertexEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::FlipVertexEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::SurfaceEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::PolyhedronSplitEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::SplitMergeEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::EdgeSplitEvent>(event)->getNode();
        } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
            node = std::dynamic_pointer_cast<data::_3d::skel::PierceEvent>(event)->getNode();
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
                result = data::_3d::skel::ConstOffsetEvent::create();
            } else if (etype == AbstractEvent::EDGE_EVENT) {
                data::_3d::skel::EdgeEventSPtr event =
                        data::_3d::skel::EdgeEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::EDGE_MERGE_EVENT) {
                data::_3d::skel::EdgeMergeEventSPtr event =
                        data::_3d::skel::EdgeMergeEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::TRIANGLE_EVENT) {
                data::_3d::skel::TriangleEventSPtr event =
                        data::_3d::skel::TriangleEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
                data::_3d::skel::DblEdgeMergeEventSPtr event =
                        data::_3d::skel::DblEdgeMergeEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::DBL_TRIANGLE_EVENT) {
                data::_3d::skel::DblTriangleEventSPtr event =
                        data::_3d::skel::DblTriangleEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::TETRAHEDRON_EVENT) {
                data::_3d::skel::TetrahedronEventSPtr event =
                        data::_3d::skel::TetrahedronEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::VERTEX_EVENT) {
                data::_3d::skel::VertexEventSPtr event =
                        data::_3d::skel::VertexEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::FLIP_VERTEX_EVENT) {
                data::_3d::skel::FlipVertexEventSPtr event =
                        data::_3d::skel::FlipVertexEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::SURFACE_EVENT) {
                data::_3d::skel::SurfaceEventSPtr event =
                        data::_3d::skel::SurfaceEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
                data::_3d::skel::PolyhedronSplitEventSPtr event =
                        data::_3d::skel::PolyhedronSplitEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::SPLIT_MERGE_EVENT) {
                data::_3d::skel::SplitMergeEventSPtr event =
                        data::_3d::skel::SplitMergeEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::EDGE_SPLIT_EVENT) {
                data::_3d::skel::EdgeSplitEventSPtr event =
                        data::_3d::skel::EdgeSplitEvent::create();
                event->setNode(node);
                result = event;
            } else if (etype == AbstractEvent::PIERCE_EVENT) {
                data::_3d::skel::PierceEventSPtr event =
                        data::_3d::skel::PierceEvent::create();
                event->setNode(node);
                result = event;
            } else {
                std::cout << "Error: etype=" << etype
                          << " does not exist." << std::endl;
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
