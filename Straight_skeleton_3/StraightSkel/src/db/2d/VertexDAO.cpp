/**
 * @file   db/2d/VertexDAO.cpp
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#include "db/2d/VertexDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/2d/PointDAO.h"
#include "db/2d/PolygonDAO.h"

namespace db { namespace _2d {

VertexDAO::VertexDAO() {
    // intentionally does nothing
}

VertexDAO::~VertexDAO() {
    // intentionally does nothing
}

std::string VertexDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Vertices (\n"
            "  PolyID INTEGER NOT NULL,\n"
            "  VID INTEGER NOT NULL,\n"
            "  PointID INTEGER,\n"
            "  PRIMARY KEY (PolyID, VID)\n"
            ");");
    return schema;
}

int VertexDAO::nextVID(int polyid) {
    int vid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(VID) FROM Vertices WHERE PolyID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        vid = 1;
        if (stmt->execute()) {
            vid = stmt->getInteger(0) + 1;
        }
    }
    return vid;
}

int VertexDAO::insert(VertexSPtr vertex) {
    int result = -1;
    if (!vertex->getPolygon()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyid = vertex->getPolygon()->getID();
    if (polyid <= 0) {
        PolygonDAOSPtr dao_polygon = DAOFactory::getPolygonDAO();
        polyid = dao_polygon->createPolyID(vertex->getPolygon());
    }
    if (polyid > 0) {
        PointDAOSPtr dao_point = DAOFactory::getPointDAO();
        int pointid = dao_point->insert(vertex->getPoint());
        if (pointid > 0) {
            int vid = nextVID(polyid);
            std::string sql("INSERT INTO Vertices (PolyID, VID, PointID) "
                    "VALUES (?, ?, ?);");
            SQLiteStmtSPtr stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyid);
                stmt->bindInteger(2, vid);
                stmt->bindInteger(3, pointid);
                if (stmt->execute() > 0) {
                    vertex->setID(vid);
                    result = vid;
                }
            }
        }
    }
    return result;
}

bool VertexDAO::del(VertexSPtr vertex) {
    bool result = false;
    if (!vertex->getPolygon()) {
        return false;
    }
    int polyid = vertex->getPolygon()->getID();
    int vid = vertex->getID();
    if (polyid < 0 || vid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Vertices WHERE PolyID=? AND VID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        stmt->bindInteger(2, vid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    sql = "DELETE FROM Edges WHERE PolyID=? AND (VID_SRC=? OR VID_DST=?);";
    stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        stmt->bindInteger(2, vid);
        stmt->bindInteger(3, vid);
        stmt->execute();
    }
    if (result) {
        vertex->setID(-1);
    }
    return result;
}

VertexSPtr VertexDAO::find(int polyid, int vid) {
    VertexSPtr result = VertexSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT PointID FROM Vertices WHERE PolyID=? AND VID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        stmt->bindInteger(2, vid);
        if (stmt->execute() > 0) {
            int point_id = stmt->getInteger(0);
            PointDAOSPtr dao_point = DAOFactory::getPointDAO();
            Point2SPtr point = dao_point->find(point_id);
            result = Vertex::create(point);
            result->setID(vid);
        }
    }
    return result;
}

bool VertexDAO::update(VertexSPtr vertex) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
