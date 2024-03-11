/**
 * @file   db/3d/NodeDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-03
 */

#include "db/3d/NodeDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/PointDAO.h"
#include "db/3d/StraightSkeletonDAO.h"

namespace db { namespace _3d {

NodeDAO::NodeDAO() {
    // intentionally does nothing
}

NodeDAO::~NodeDAO() {
    // intentionally does nothing
}

std::string NodeDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Nodes (\n"
            "  SkelID INTEGER NOT NULL,\n"
            "  NID INTEGER NOT NULL,\n"
            "  PointID INTEGER,\n"
            "  offset REAL,\n"
            "  PRIMARY KEY (SkelID, NID)\n"
            ");");
    return schema;
}

int NodeDAO::nextNID(int skelid) {
    int nid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(NID) FROM Nodes WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        nid = 1;
        if (stmt->execute()) {
            nid = stmt->getInteger(0) + 1;
        }
    }
    return nid;
}

int NodeDAO::insert(NodeSPtr node) {
    int result = -1;
    if (!node->getSkel()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int skelid = node->getSkel()->getID();
    if (skelid <= 0) {
        StraightSkeletonDAOSPtr dao_skel = DAOFactory::getStraightSkeletonDAO();
        skelid = dao_skel->createSkelID(node->getSkel());
    }
    if (skelid > 0) {
        PointDAOSPtr dao_point = DAOFactory::getPointDAO();
        int pointid = dao_point->insert(node->getPoint());
        if (pointid > 0) {
            int nid = nextNID(skelid);
            std::string sql("INSERT INTO Nodes (SkelID, NID, PointID, offset) "
                    "VALUES (?, ?, ?, ?);");
            SQLiteStmtSPtr stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->bindInteger(2, nid);
                stmt->bindInteger(3, pointid);
                stmt->bindDouble(4, node->getOffset());
                if (stmt->execute() > 0) {
                    node->setID(nid);
                    result = nid;
                }
            }
        }
    }
    return result;
}

bool NodeDAO::del(NodeSPtr node) {
    bool result = false;
    if (!node->getSkel()) {
        return false;
    }
    int skelid = node->getSkel()->getID();
    int nid = node->getID();
    if (skelid < 0 || nid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Nodes WHERE SkelID=? AND NID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, nid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    sql = "DELETE FROM Arcs WHERE SkelID=? AND (NID_SRC=? OR NID_DST=?);";
    stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, nid);
        stmt->bindInteger(3, nid);
        stmt->execute();
    }
    if (result) {
        node->setID(-1);
    }
    return result;
}

NodeSPtr NodeDAO::find(int skelid, int nid) {
    NodeSPtr result = NodeSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT PointID, offset FROM Nodes WHERE SkelID=? AND NID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, nid);
        if (stmt->execute() > 0) {
            int point_id = stmt->getInteger(0);
            double offset = stmt->getDouble(1);
            PointDAOSPtr dao_point = DAOFactory::getPointDAO();
            Point3SPtr point = dao_point->find(point_id);
            result = Node::create(point);
            result->setOffset(offset);
            result->setID(nid);
        }
    }
    return result;
}

bool NodeDAO::update(NodeSPtr node) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
