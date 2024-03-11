/**
 * @file   db/3d/ArcDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-03
 */

#include "db/3d/ArcDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/NodeDAO.h"
#include "db/3d/StraightSkeletonDAO.h"

namespace db { namespace _3d {

ArcDAO::ArcDAO() {
    // intentionally does nothing
}

ArcDAO::~ArcDAO() {
    // intentionally does nothing
}

std::string ArcDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Arcs (\n"
            "  SkelID INTEGER NOT NULL,\n"
            "  AID INTEGER NOT NULL,\n"
            "  NID_SRC INTEGER,\n"
            "  NID_DST INTEGER,\n"
            "  PRIMARY KEY (SkelID, AID)\n"
            ");");
    return schema;
}

int ArcDAO::nextAID(int skelid) {
    int aid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(AID) FROM Arcs WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        aid = 1;
        if (stmt->execute()) {
            aid = stmt->getInteger(0) + 1;
        }
    }
    return aid;
}

int ArcDAO::insert(ArcSPtr arc) {
    int result = -1;
    if (!arc->getSkel()) {
        return -1;
    }
    if (!arc->hasNodeDst()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int skelid = arc->getSkel()->getID();
    if (skelid <= 0) {
        StraightSkeletonDAOSPtr dao_skel = DAOFactory::getStraightSkeletonDAO();
        skelid = dao_skel->createSkelID(arc->getSkel());
    }
    if (skelid > 0) {
        int aid = nextAID(skelid);
        NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
        if (arc->getNodeSrc()->getID() > 0) {
            dao_node->update(arc->getNodeSrc());
        } else {
            dao_node->insert(arc->getNodeSrc());
        }
        if (arc->getNodeDst()->getID() > 0) {
            dao_node->update(arc->getNodeDst());
        } else {
            dao_node->insert(arc->getNodeDst());
        }
        std::string sql("INSERT INTO Arcs (SkelID, AID, NID_SRC, NID_DST) "
                "VALUES (?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, skelid);
            stmt->bindInteger(2, aid);
            stmt->bindInteger(3, arc->getNodeSrc()->getID());
            stmt->bindInteger(4, arc->getNodeDst()->getID());
            if (stmt->execute() > 0) {
                arc->setID(aid);
                result = aid;
            }
        }
    }
    return result;
}

bool ArcDAO::del(ArcSPtr arc) {
    bool result = false;
    if (!arc->getSkel()) {
        return false;
    }
    int skelid = arc->getSkel()->getID();
    int aid = arc->getID();
    if (skelid < 0 || aid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Arcs WHERE SkelID=? AND AID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, aid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    if (result) {
        arc->setID(-1);
    }
    return result;
}

ArcSPtr ArcDAO::find(int skelid, int aid) {
    ArcSPtr result = ArcSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT NID_SRC, NID_DST FROM Arcs WHERE SkelID=? AND AID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, aid);
        if (stmt->execute() > 0) {
            int nid_src = stmt->getInteger(0);
            int nid_dst = stmt->getInteger(1);
            NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
            NodeSPtr node_src = dao_node->find(skelid, nid_src);
            NodeSPtr node_dst = dao_node->find(skelid, nid_dst);
            result = Arc::create(node_src, node_dst);
            result->setID(aid);
            node_src->addArc(result);
            node_dst->addArc(result);
        }
    }
    return result;
}

bool ArcDAO::update(ArcSPtr arc) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
