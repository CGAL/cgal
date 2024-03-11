/**
 * @file   db/3d/StraightSkeletonDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-03
 */

#include "db/3d/StraightSkeletonDAO.h"

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
#include "db/3d/EventDAO.h"
#include "db/3d/NodeDAO.h"
#include "db/3d/ArcDAO.h"
#include "db/3d/PolyhedronDAO.h"
#include "db/3d/SheetDAO.h"
#include <list>
#include <map>

namespace db { namespace _3d {

StraightSkeletonDAO::StraightSkeletonDAO() {
    // intentionally does nothing
}

StraightSkeletonDAO::~StraightSkeletonDAO() {
    // intentionally does nothing
}

std::string StraightSkeletonDAO::getTableSchema() const {
    std::string schema("CREATE TABLE StraightSkeletons (\n"
            "  SkelID INTEGER PRIMARY KEY,\n"
            "  PolyhedronID INTEGER,\n"
            "  config TEXT,\n"
            "  description TEXT,\n"
            "  created INTEGER\n"
            ");");
    return schema;
}

int StraightSkeletonDAO::nextSkelID() {
    int skelid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(SkelID) FROM StraightSkeletons;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        skelid = 1;
        if (stmt->execute()) {
            skelid = stmt->getInteger(0) + 1;
        }
    }
    return skelid;
}

int StraightSkeletonDAO::createSkelID(StraightSkeletonSPtr skel) {
    int result = -1;
    int skelid = nextSkelID();
    int polyhedronid = 0;
    if (skel->getPolyhedron()) {
        polyhedronid = skel->getPolyhedron()->getID();
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql;
    if (polyhedronid > 0) {
        sql = "INSERT INTO StraightSkeletons (SkelID, PolyhedronID, config, description, created) "
            "VALUES (?, ?, ?, ?, strftime('%s','now'));";
    } else {
        sql = "INSERT INTO StraightSkeletons (SkelID, config, description, created) "
            "VALUES (?, ?, ?, strftime('%s','now'));";
    }
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        if (skelid > 0) {
            stmt->bindInteger(1, skelid);
            if (polyhedronid > 0) {
                stmt->bindInteger(2, polyhedronid);
                stmt->bindString(3, skel->getConfig());
                stmt->bindString(4, skel->getDescription());
            } else {
                stmt->bindString(2, skel->getConfig());
                stmt->bindString(3, skel->getDescription());
            }
            if (stmt->execute()) {
                result = skelid;
                skel->setID(result);
            }
        }
    }
    return result;
}

int StraightSkeletonDAO::insert(StraightSkeletonSPtr skel) {
    int result = -1;
    skel->resetAllIDs();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    int skelid = createSkelID(skel);
    if (skelid > 0) {
        NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
        std::list<NodeSPtr>::iterator it_n = skel->nodes().begin();
        while (it_n != skel->nodes().end()) {
            NodeSPtr node = *it_n++;
            if (node->getID() > 0) {
                dao_node->update(node);
            } else {
                dao_node->insert(node);
            }
        }
        ArcDAOSPtr dao_arc = DAOFactory::getArcDAO();
        std::list<ArcSPtr>::iterator it_a = skel->arcs().begin();
        while (it_a != skel->arcs().end()) {
            ArcSPtr arc = *it_a++;
            if (arc->getID() > 0) {
                dao_arc->update(arc);
            } else {
                dao_arc->insert(arc);
            }
        }
        SheetDAOSPtr dao_sheet = DAOFactory::getSheetDAO();
        std::list<SheetSPtr>::iterator it_s = skel->sheets().begin();
        while (it_s != skel->sheets().end()) {
            SheetSPtr sheet = *it_s++;
            if (sheet->getID() > 0) {
                dao_sheet->update(sheet);
            } else {
                dao_sheet->insert(sheet);
            }
        }
        EventDAOSPtr dao_event = DAOFactory::getEventDAO();
        std::list<AbstractEventSPtr>::iterator it_e = skel->events().begin();
        while (it_e != skel->events().end()) {
            AbstractEventSPtr event = *it_e++;
            if (event->getID() > 0) {
                dao_event->update(event);
            } else {
                dao_event->insert(event);
            }
        }
        result = skelid;
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

bool StraightSkeletonDAO::del(StraightSkeletonSPtr skel) {
    bool result = false;
    int skelid = skel->getID();
    if (skelid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    std::string sql("DELETE FROM StraightSkeletons WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        if (stmt->execute() > 0) {
            sql = "DELETE FROM Events WHERE SkelID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->execute();
            }
            sql = "DELETE FROM Sheets WHERE SkelID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->execute();
            }
            sql = "DELETE FROM Arcs WHERE SkelID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->execute();
            }
            sql = "DELETE FROM Nodes WHERE SkelID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->execute();
            }
            result = true;
            skel->setID(-1);
        }
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

StraightSkeletonSPtr StraightSkeletonDAO::find(int skelid) {
    StraightSkeletonSPtr result = StraightSkeletonSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
    EventDAOSPtr dao_event = DAOFactory::getEventDAO();
    std::string sql("SELECT SkelID, config FROM StraightSkeletons WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        if (stmt->execute() > 0) {
            result = StraightSkeleton::create();
            result->setID(skelid);
            std::string config = stmt->getString(1);
            result->setConfig(config);
            std::map<int, NodeSPtr> nodes;
            sql = "SELECT NID FROM Nodes WHERE SkelID=? ORDER BY NID ASC;";
            SQLiteStmtSPtr stmt_n = db->prepare(sql);
            if (stmt_n) {
                stmt_n->bindInteger(1, skelid);
                stmt_n->execute();
                while (stmt_n->fetchRow()) {
                    int nid = stmt_n->getInteger(0);
                    NodeSPtr node = dao_node->find(skelid, nid);
                    nodes[nid] = node;
                    result->addNode(node);
                }
                stmt_n->close();
            }
            std::map<int, ArcSPtr> arcs;
            sql = "SELECT AID, NID_SRC, NID_DST FROM Arcs WHERE SkelID=? ORDER BY AID ASC;";
            SQLiteStmtSPtr stmt_a = db->prepare(sql);
            if (stmt_a) {
                stmt_a->bindInteger(1, skelid);
                stmt_a->execute();
                while (stmt_a->fetchRow()) {
                    int aid = stmt_a->getInteger(0);
                    int nid_src = stmt_a->getInteger(1);
                    int nid_dst = stmt_a->getInteger(2);
                    ArcSPtr arc = Arc::create(nodes[nid_src], nodes[nid_dst]);
                    arc->setID(aid);
                    nodes[nid_src]->addArc(arc);
                    nodes[nid_dst]->addArc(arc);
                    arcs[aid] = arc;
                    result->addArc(arc);
                }
                stmt_a->close();
            }
            sql = "SELECT SID FROM Sheets WHERE SkelID=? ORDER BY SID ASC;";
            SQLiteStmtSPtr stmt_s = db->prepare(sql);
            if (stmt_s) {
                stmt_s->bindInteger(1, skelid);
                stmt_s->execute();
                while (stmt_s->fetchRow()) {
                    int sid = stmt_s->getInteger(0);
                    SheetSPtr sheet = Sheet::create();
                    sheet->setID(sid);
                    std::string sql_sa = "SELECT AID FROM Sheets_Arcs "
                            "WHERE SkelID=? AND SID=? ORDER BY AID ASC;";
                    SQLiteStmtSPtr stmt_sa = db->prepare(sql_sa);
                    if (stmt_sa) {
                        stmt_sa->bindInteger(1, skelid);
                        stmt_sa->bindInteger(2, sid);
                        stmt_sa->execute();
                        while (stmt_sa->fetchRow()) {
                            int aid = stmt_sa->getInteger(0);
                            ArcSPtr arc = arcs[aid];
                            sheet->addNode(arc->getNodeSrc());
                            sheet->addArc(arc);
                        }
                        stmt_sa->close();
                    }
                    result->addSheet(sheet);
                }
                stmt_s->close();
            }
            sql = "SELECT EventID, NID FROM Events WHERE SkelID=? ORDER BY EventID ASC;";
            SQLiteStmtSPtr stmt_e = db->prepare(sql);
            if (stmt_e) {
                stmt_e->bindInteger(1, skelid);
                stmt_e->execute();
                while (stmt_e->fetchRow()) {
                    int eventid = stmt_e->getInteger(0);
                    int nid = stmt_e->getInteger(1);
                    AbstractEventSPtr event = dao_event->find(skelid, eventid);
                    if (event->getType() == AbstractEvent::EDGE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::EdgeEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::EdgeMergeEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::TriangleEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::DblEdgeMergeEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::DblTriangleEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::TetrahedronEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::VertexEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::FlipVertexEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::SurfaceEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::PolyhedronSplitEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::SplitMergeEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::EdgeSplitEvent>(event)->setNode(nodes[nid]);
                    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
                        std::dynamic_pointer_cast<data::_3d::skel::PierceEvent>(event)->setNode(nodes[nid]);
                    }
                    result->addEvent(event);
                }
                stmt_e->close();
            }
        }
    }
    return result;
}

int StraightSkeletonDAO::findPolyhedronID(int skelid) {
    int result = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT PolyhedronID FROM StraightSkeletons WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        if (stmt->execute() > 0) {
            result = stmt->getInteger(0);
        }
    }
    return result;
}

bool StraightSkeletonDAO::update(StraightSkeletonSPtr skel) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
