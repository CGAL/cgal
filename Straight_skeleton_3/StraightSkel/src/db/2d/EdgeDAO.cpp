/**
 * @file   db/2d/EdgeDAO.cpp
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#include "db/2d/EdgeDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/2d/VertexDAO.h"
#include "db/2d/PolygonDAO.h"

namespace db { namespace _2d {

EdgeDAO::EdgeDAO() {
    // intentionally does nothing
}

EdgeDAO::~EdgeDAO() {
    // intentionally does nothing
}

std::string EdgeDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Edges (\n"
            "  PolyID INTEGER NOT NULL,\n"
            "  EID INTEGER NOT NULL,\n"
            "  VID_SRC INTEGER,\n"
            "  VID_DST INTEGER,\n"
            "  PRIMARY KEY (PolyID, EID)\n"
            ");");
    return schema;
}

std::string EdgeDAO::getTable2Schema() const {
    std::string schema("CREATE TABLE SkelEdgeData (\n"
            "  PolyID INTEGER NOT NULL,\n"
            "  EID INTEGER NOT NULL,\n"
            "  speed REAL,\n"
            "  PRIMARY KEY (PolyID, EID)\n"
            ");");
    return schema;
}

int EdgeDAO::nextEID(int polyid) {
    int eid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(EID) FROM Edges WHERE PolyID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        eid = 1;
        if (stmt->execute()) {
            eid = stmt->getInteger(0) + 1;
        }
    }
    return eid;
}

int EdgeDAO::insert(EdgeSPtr edge) {
    int result = -1;
    if (!edge->getPolygon()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyid = edge->getPolygon()->getID();
    if (polyid <= 0) {
        PolygonDAOSPtr dao_polygon = DAOFactory::getPolygonDAO();
        polyid = dao_polygon->createPolyID(edge->getPolygon());
    }
    if (polyid > 0) {
        int eid = nextEID(polyid);
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
        std::string sql("INSERT INTO Edges (PolyID, EID, VID_SRC, VID_DST) "
                "VALUES (?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, polyid);
            stmt->bindInteger(2, eid);
            stmt->bindInteger(3, edge->getVertexSrc()->getID());
            stmt->bindInteger(4, edge->getVertexDst()->getID());
            if (stmt->execute() > 0) {
                edge->setID(eid);
                result = eid;
                if (edge->hasData()) {
                    SkelEdgeDataSPtr data = std::dynamic_pointer_cast<SkelEdgeData>(
                            edge->getData());
                    sql = "INSERT INTO SkelEdgeData (PolyID, EID, speed) "
                            "VALUES (?, ?, ?);";
                    stmt = db->prepare(sql);
                    if (stmt) {
                        stmt->bindInteger(1, polyid);
                        stmt->bindInteger(2, eid);
                        stmt->bindDouble(3, data->getSpeed());
                        stmt->execute();
                    }
                }
            }
        }
    }
    return result;
}

bool EdgeDAO::del(EdgeSPtr edge) {
    bool result = false;
    if (!edge->getPolygon()) {
        return false;
    }
    int polyid = edge->getPolygon()->getID();
    int eid = edge->getID();
    if (polyid < 0 || eid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Edges WHERE PolyID=? AND EID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        stmt->bindInteger(2, eid);
        if (stmt->execute() > 0) {
            result = true;
            if (edge->hasData()) {
                sql = "DELETE FROM SkelEdgeData WHERE PolyID=? AND EID=?;";
                stmt = db->prepare(sql);
                if (stmt) {
                    stmt->bindInteger(1, polyid);
                    stmt->bindInteger(2, eid);
                    stmt->execute();
                }
            }

        }
    }
    if (result) {
        edge->setID(-1);
    }
    return result;
}

EdgeSPtr EdgeDAO::find(int polyid, int eid) {
    EdgeSPtr result = EdgeSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT VID_SRC, VID_DST FROM Edges WHERE PolyID=? AND EID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        stmt->bindInteger(2, eid);
        if (stmt->execute() > 0) {
            int vid_src = stmt->getInteger(0);
            int vid_dst = stmt->getInteger(1);
            VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
            VertexSPtr vertex_src = dao_vertex->find(polyid, vid_src);
            VertexSPtr vertex_dst = dao_vertex->find(polyid, vid_dst);
            result = Edge::create(vertex_src, vertex_dst);
            result->setID(eid);
            sql = "SELECT speed FROM SkelEdgeData WHERE PolyID=? AND EID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyid);
                stmt->bindInteger(2, eid);
                if (stmt->execute() > 0) {
                    double speed = stmt->getDouble(0);
                    if (speed != 0.0 && speed != 1.0) {
                        SkelEdgeDataSPtr data = SkelEdgeData::create(result);
                        data->setSpeed(speed);
                    }
                }
            }
        }
    }
    return result;
}

bool EdgeDAO::update(EdgeSPtr edge) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
