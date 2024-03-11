/**.
 * @file   db/3d/SheetDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-05
 */

#include "db/3d/SheetDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/NodeDAO.h"
#include "db/3d/ArcDAO.h"
#include "db/3d/StraightSkeletonDAO.h"
#include <list>
#include <map>

namespace db { namespace _3d {

SheetDAO::SheetDAO() {
    // intentionally does nothing
}

SheetDAO::~SheetDAO() {
    // intentionally does nothing
}

std::string SheetDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Sheets (\n"
            "  SkelID INTEGER NOT NULL,\n"
            "  SID INTEGER NOT NULL,\n"
            "  PRIMARY KEY (SkelID, SID)\n"
            ");");
    return schema;
}

std::string SheetDAO::getTable2Schema() const {
    std::string schema("CREATE TABLE Sheets_Arcs (\n"
            "  SkelID INTEGER NOT NULL,\n"
            "  SID INTEGER NOT NULL,\n"
            "  AID INTEGER NOT NULL,\n"
            "  PRIMARY KEY (SkelID, SID, AID)\n"
            ");");
    return schema;
}

int SheetDAO::nextSID(int skelid) {
    int sid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(SID) FROM Sheets WHERE SkelID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        sid = 1;
        if (stmt->execute()) {
            sid = stmt->getInteger(0) + 1;
        }
    }
    return sid;
}

int SheetDAO::insert(SheetSPtr sheet) {
    int result = -1;
    if (!sheet->getSkel()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int skelid = sheet->getSkel()->getID();
    if (skelid <= 0) {
        StraightSkeletonDAOSPtr dao_skel = DAOFactory::getStraightSkeletonDAO();
        skelid = dao_skel->createSkelID(sheet->getSkel());
    }
    if (skelid > 0) {
        int sid = nextSID(skelid);
        std::string sql("INSERT INTO Sheets (SkelID, SID) VALUES (?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, skelid);
            stmt->bindInteger(2, sid);
            if (stmt->execute() > 0) {
                sheet->setID(sid);
                result = sid;

                ArcDAOSPtr dao_arc = DAOFactory::getArcDAO();
                sql = "INSERT INTO Sheets_Arcs (SkelID, SID, AID) "
                        "VALUES (?, ?, ?);";
                SQLiteStmtSPtr stmt_sa = db->prepare(sql);
                if (stmt_sa) {
                    std::list<ArcSPtr>::iterator it_a = sheet->arcs().begin();
                    while (it_a != sheet->arcs().end()) {
                        ArcSPtr arc = *it_a++;
                        if (arc->getID() > 0) {
                            dao_arc->update(arc);
                        } else {
                            dao_arc->insert(arc);
                        }
                        stmt_sa->bindInteger(1, skelid);
                        stmt_sa->bindInteger(2, sid);
                        stmt_sa->bindInteger(3, arc->getID());
                        stmt_sa->execute();
                        stmt_sa->reset();
                    }
                    stmt_sa->close();
                }
            }
        }
    }
    return result;
}

bool SheetDAO::del(SheetSPtr sheet) {
    int result = false;
    if (!sheet->getSkel()) {
        return false;
    }
    int skelid = sheet->getSkel()->getID();
    int sid = sheet->getID();
    if (skelid < 0 || sid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Sheets WHERE SkelID=? AND SID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, sid);
        if (stmt->execute() > 0) {
            sql = "DELETE FROM Sheets_Arcs WHERE SkelID=? AND SID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, skelid);
                stmt->bindInteger(2, sid);
                stmt->execute();
            }
            result = true;
            sheet->setID(-1);
        }
    }
    return result;
}

SheetSPtr SheetDAO::find(int skelid, int sid) {
    SheetSPtr result = SheetSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    NodeDAOSPtr dao_node = DAOFactory::getNodeDAO();
    std::string sql("SELECT SID FROM Sheets "
            "WHERE SkelID=? AND SID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, skelid);
        stmt->bindInteger(2, sid);
        if (stmt->execute() > 0) {
            result = Sheet::create();
            result->setID(sid);
            std::map<int, NodeSPtr> nodes;
            sql = "SELECT Arcs.AID, NID_SRC, NID_DST FROM Arcs "
                    "JOIN Sheets_Arcs ON (Arcs.SkelID=Sheets_Arcs.SkelID AND Arcs.AID=Sheets_Arcs.AID) "
                    "WHERE Arcs.SkelID=? AND SID=?";
            SQLiteStmtSPtr stmt_a = db->prepare(sql);
            if (stmt_a) {
                stmt_a->bindInteger(1, skelid);
                stmt_a->bindInteger(2, sid);
                stmt_a->execute();
                while (stmt_a->fetchRow()) {
                    int aid = stmt_a->getInteger(0);
                    int nid_src = stmt_a->getInteger(1);
                    int nid_dst = stmt_a->getInteger(2);
                    if (!nodes[nid_src]) {
                        NodeSPtr node = dao_node->find(skelid, nid_src);
                        nodes[nid_src] = node;
                        result->addNode(node);
                    }
                    if (!nodes[nid_dst]) {
                        NodeSPtr node = dao_node->find(skelid, nid_dst);
                        nodes[nid_dst] = node;
                        result->addNode(node);
                    }
                    ArcSPtr arc = Arc::create(nodes[nid_src], nodes[nid_dst]);
                    arc->setID(aid);
                    result->addArc(arc);
                }
                stmt_a->close();
            }
        }
    }
    return result;
}

bool SheetDAO::update(SheetSPtr sheet) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
