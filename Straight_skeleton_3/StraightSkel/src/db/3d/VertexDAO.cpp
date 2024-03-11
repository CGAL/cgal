#include "db/3d/VertexDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/PointDAO.h"
#include "db/3d/PolyhedronDAO.h"

namespace db { namespace _3d {

VertexDAO::VertexDAO() {
    // intentionally does nothing
}

VertexDAO::~VertexDAO() {
    // intentionally does nothing
}

std::string VertexDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Vertices (\n"
            "  PolyhedronID INTEGER NOT NULL,\n"
            "  VID INTEGER NOT NULL,\n"
            "  PointID INTEGER,\n"
            "  PRIMARY KEY (PolyhedronID, VID)\n"
            ");");
    return schema;
}

int VertexDAO::nextVID(int polyhedronid) {
    int vid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(VID) FROM Vertices WHERE PolyhedronID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        vid = 1;
        if (stmt->execute()) {
            vid = stmt->getInteger(0) + 1;
        }
    }
    return vid;
}

int VertexDAO::insert(VertexSPtr vertex) {
    int result = -1;
    if (!vertex->getPolyhedron()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyhedronid = vertex->getPolyhedron()->getID();
    if (polyhedronid <= 0) {
        PolyhedronDAOSPtr dao_polyhedron = DAOFactory::getPolyhedronDAO();
        polyhedronid = dao_polyhedron->createPolyhedronID(vertex->getPolyhedron());
    }
    if (polyhedronid > 0) {
        PointDAOSPtr dao_point = DAOFactory::getPointDAO();
        int pointid = dao_point->insert(vertex->getPoint());
        if (pointid > 0) {
            int vid = nextVID(polyhedronid);
            std::string sql("INSERT INTO Vertices (PolyhedronID, VID, PointID) "
                    "VALUES (?, ?, ?);");
            SQLiteStmtSPtr stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedronid);
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
    if (!vertex->getPolyhedron()) {
        return false;
    }
    int polyhedronid = vertex->getPolyhedron()->getID();
    int vid = vertex->getID();
    if (polyhedronid < 0 || vid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Vertices WHERE PolyhedronID=? AND VID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, vid);
        if (stmt->execute() > 0) {
            result = true;
        }
    }
    sql = "DELETE FROM Edges WHERE PolyhedronID=? AND "
            "(VID_SRC=? OR VID_DST=?);";
    stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, vid);
        stmt->bindInteger(3, vid);
        stmt->execute();
    }
    sql = "DELETE FROM Triangles WHERE PolyhedronID=? AND "
            "(VID_1=? OR VID_2=? OR VID_3=?);";
    stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, vid);
        stmt->bindInteger(3, vid);
        stmt->bindInteger(4, vid);
        stmt->execute();
    }
    if (result) {
        vertex->setID(-1);
    }
    return result;
}

VertexSPtr VertexDAO::find(int polyhedronid, int vid) {
    VertexSPtr result = VertexSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT PointID FROM Vertices WHERE PolyhedronID=? AND VID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, vid);
        if (stmt->execute() > 0) {
            int point_id = stmt->getInteger(0);
            PointDAOSPtr dao_point = DAOFactory::getPointDAO();
            Point3SPtr point = dao_point->find(point_id);
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
