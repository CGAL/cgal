#include "db/3d/TriangleDAO.h"

#include "data/3d/Triangle.h"
#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/VertexDAO.h"
#include "db/3d/FacetDAO.h"
#include "db/3d/PolyhedronDAO.h"

namespace db { namespace _3d {

TriangleDAO::TriangleDAO() {
    // intentionally does nothing
}

TriangleDAO::~TriangleDAO() {
    // intentionally does nothing
}

std::string TriangleDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Triangles (\n"
            "  PolyhedronID INTEGER NOT NULL,\n"
            "  FID INTEGER NOT NULL,\n"
            "  TID INTEGER NOT NULL,\n"
            "  VID_1 INTEGER,\n"
            "  VID_2 INTEGER,\n"
            "  VID_3 INTEGER,\n"
            "  PRIMARY KEY (PolyhedronID, FID, TID)\n"
            ");");
    return schema;
}

int TriangleDAO::nextTID(int polyhedronid, int fid) {
    int tid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(TID) FROM Triangles "
            "WHERE PolyhedronID=? AND FID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, fid);
        tid = 1;
        if (stmt->execute()) {
            tid = stmt->getInteger(0) + 1;
        }
    }
    return tid;
}

int TriangleDAO::insert(TriangleSPtr triangle) {
    int result = -1;
    if (!triangle->getFacet()) {
        return -1;
    }
    if (!triangle->getFacet()->getPolyhedron()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyhedronid = triangle->getFacet()->getPolyhedron()->getID();
    if (polyhedronid <= 0) {
        PolyhedronDAOSPtr dao_polyhedron = DAOFactory::getPolyhedronDAO();
        polyhedronid = dao_polyhedron->createPolyhedronID(
                triangle->getFacet()->getPolyhedron());
    }
    int fid = triangle->getFacet()->getID();
    if (fid <= 0) {
        FacetDAOSPtr dao_facet = DAOFactory::getFacetDAO();
        fid = dao_facet->createFID(triangle->getFacet());
    }
    if (polyhedronid > 0 && fid > 0) {
        int tid = nextTID(polyhedronid, fid);
        VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
        VertexSPtr verts[3];
        for (unsigned int i = 0; i < 3; i++) {
            VertexSPtr vertex = triangle->getVertex(i);
            verts[i] = vertex;
            if (vertex->getID() > 0) {
                dao_vertex->update(vertex);
            } else {
                dao_vertex->insert(vertex);
            }
        }
        std::string sql("INSERT INTO Triangles (PolyhedronID, FID, TID, VID_1, VID_2, VID_3) "
            "VALUES (?, ?, ?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, polyhedronid);
            stmt->bindInteger(2, fid);
            stmt->bindInteger(3, tid);
            stmt->bindInteger(4, verts[0]->getID());
            stmt->bindInteger(5, verts[1]->getID());
            stmt->bindInteger(6, verts[2]->getID());
            if (stmt->execute() > 0) {
                triangle->setID(tid);
                result = tid;
            }
        }
    }
    return result;
}

bool TriangleDAO::del(TriangleSPtr triangle) {
    bool result = false;
    if (!triangle->getFacet()) {
        return false;
    }
    if (!triangle->getFacet()->getPolyhedron()) {
        return false;
    }
    int tid = triangle->getID();
    int fid = triangle->getFacet()->getID();
    int polyhedronid = triangle->getFacet()->getPolyhedron()->getID();
    if (polyhedronid < 0 || fid < 0 || tid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Triangles WHERE PolyhedronID=? AND FID=? AND TID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, fid);
        stmt->bindInteger(3, tid);
        if (stmt->execute() > 0) {
            triangle->setID(-1);
            result = true;
        }
    }
    return result;
}

TriangleSPtr TriangleDAO::find(int polyhedronid, int fid, int tid) {
    TriangleSPtr result = TriangleSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT VID_1, VID_2, VID_3 FROM Triangles "
            "WHERE PolyhedronID=? AND FID=? AND TID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, fid);
        stmt->bindInteger(3, tid);
        if (stmt->execute() > 0) {
            int vids[3];
            for (unsigned int i = 0; i < 3; i++) {
                vids[i] = stmt->getInteger(i);
            }
            VertexDAOSPtr dao_vertex= DAOFactory::getVertexDAO();
            VertexSPtr verts[3];
            for (unsigned int i = 0; i < 3; i++) {
                verts[i] = dao_vertex->find(polyhedronid, vids[i]);
            }
            FacetSPtr facet = Facet::create();
            facet->setID(fid);
            result = Triangle::create(facet, verts);
            result->setID(tid);
        }
    }
    return result;
}

bool TriangleDAO::update(TriangleSPtr triangle) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
