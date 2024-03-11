#include "db/3d/EdgeDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/PointDAO.h"
#include "db/3d/VertexDAO.h"
#include "db/3d/PolyhedronDAO.h"
#include "db/3d/FacetDAO.h"

namespace db { namespace _3d {

EdgeDAO::EdgeDAO() {
    // intentionally does nothing
}

EdgeDAO::~EdgeDAO() {
    // intentionally does nothing
}

std::string EdgeDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Edges (\n"
            "  PolyhedronID INTEGER NOT NULL,\n"
            "  EID INTEGER NOT NULL,\n"
            "  VID_SRC INTEGER,\n"
            "  VID_DST INTEGER,\n"
            "  FID_L INTEGER,\n"
            "  FID_R INTEGER,\n"
            "  PRIMARY KEY (PolyhedronID, EID)\n"
            ");");
    return schema;
}

int EdgeDAO::nextEID(int polyhedronid) {
    int eid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(EID) FROM Edges WHERE PolyhedronID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        eid = 1;
        if (stmt->execute()) {
            eid = stmt->getInteger(0) + 1;
        }
    }
    return eid;
}

int EdgeDAO::insert(EdgeSPtr edge) {
    int result = -1;
    if (!edge->getPolyhedron()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyhedronid = edge->getPolyhedron()->getID();
    if (polyhedronid <= 0) {
        PolyhedronDAOSPtr dao_polyhedron = DAOFactory::getPolyhedronDAO();
        polyhedronid = dao_polyhedron->createPolyhedronID(edge->getPolyhedron());
    }
    if (polyhedronid > 0) {
        int eid = nextEID(polyhedronid);
        VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
        if (edge->getVertexSrc()->getID() > 0) {
            dao_vertex->update(edge->getVertexSrc());
        } else {
            dao_vertex->insert(edge->getVertexSrc());
        }
        if (edge->getVertexDst()->getID() > 0) {
            dao_vertex->update(edge->getVertexDst());
        } else {
            dao_vertex->insert(edge->getVertexDst());
        }
        std::string sql("INSERT INTO Edges (PolyhedronID, EID, VID_SRC, VID_DST) "
                "VALUES (?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, polyhedronid);
            stmt->bindInteger(2, eid);
            stmt->bindInteger(3, edge->getVertexSrc()->getID());
            stmt->bindInteger(4, edge->getVertexDst()->getID());
            if (stmt->execute() > 0) {
                edge->setID(eid);
                result = eid;
            }
        }
    }
    return result;
}

bool EdgeDAO::del(EdgeSPtr edge) {
    bool result = false;
    if (!edge->getPolyhedron()) {
        return false;
    }
    int polyhedronid = edge->getPolyhedron()->getID();
    int eid = edge->getID();
    if (polyhedronid < 0 || eid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Edges WHERE PolyhedronID=? AND EID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, eid);
        if (stmt->execute() > 0) {
            edge->setID(-1);
            result = true;
        }
    }
    return result;
}

EdgeSPtr EdgeDAO::find(int polyhedronid, int eid) {
    EdgeSPtr result = EdgeSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT VID_SRC, VID_DST FROM Edges WHERE PolyhedronID=? AND EID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, eid);
        if (stmt->execute() > 0) {
            int vid_src = stmt->getInteger(0);
            int vid_dst = stmt->getInteger(1);
            VertexDAOSPtr dao_vertex= DAOFactory::getVertexDAO();
            VertexSPtr vertex_src = dao_vertex->find(polyhedronid, vid_src);
            VertexSPtr vertex_dst = dao_vertex->find(polyhedronid, vid_dst);
            result = Edge::create(vertex_src, vertex_dst);
            result->setID(eid);
        }
    }
    return result;
}

bool EdgeDAO::update(EdgeSPtr edge) {
    bool result = false;
    if (!edge->getPolyhedron()) {
        return false;
    }
    int polyhedronid = edge->getPolyhedron()->getID();
    int eid = edge->getID();
    if (polyhedronid < 0 || eid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("UPDATE Edges "
            "SET VID_SRC=?, VID_DST=?, FID_L=?, FID_R=? "
            "WHERE PolyhedronID=? AND EID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, edge->getVertexSrc()->getID());
        stmt->bindInteger(2, edge->getVertexDst()->getID());
        if (edge->getFacetL()) {
            if (edge->getFacetL()->getID() >= 0) {
                stmt->bindInteger(3, edge->getFacetL()->getID());
            }
        }
        if (edge->getFacetR()) {
            if (edge->getFacetR()->getID() >= 0) {
                stmt->bindInteger(4, edge->getFacetR()->getID());
            }
        }
        stmt->bindInteger(5, polyhedronid);
        stmt->bindInteger(6, eid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    return result;
}

} }
