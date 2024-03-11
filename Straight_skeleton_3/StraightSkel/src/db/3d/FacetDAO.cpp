#include "db/3d/FacetDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/VertexDAO.h"
#include "db/3d/EdgeDAO.h"
#include "db/3d/TriangleDAO.h"
#include "db/3d/PlaneDAO.h"
#include "db/3d/PolyhedronDAO.h"
#include <map>

namespace db { namespace _3d {

FacetDAO::FacetDAO() {
    // intentionally does nothing
}

FacetDAO::~FacetDAO() {
    // intentionally does nothing
}

std::string FacetDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Facets (\n"
            "  PolyhedronID INTEGER NOT NULL,\n"
            "  FID INTEGER NOT NULL,\n"
            "  PlaneID INTEGER,\n"
            "  PRIMARY KEY (PolyhedronID, FID)\n"
            ");");
    return schema;
}

int FacetDAO::nextFID(int polyhedronid) {
    int fid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(FID) FROM Facets WHERE PolyhedronID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        fid = 1;
        if (stmt->execute()) {
            fid = stmt->getInteger(0) + 1;
        }
    }
    return fid;
}

int FacetDAO::createFID(FacetSPtr facet) {
    int result = -1;
    if (!facet->getPolyhedron()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int polyhedronid = facet->getPolyhedron()->getID();
    if (polyhedronid <= 0) {
        PolyhedronDAOSPtr dao_polyhedron = DAOFactory::getPolyhedronDAO();
        polyhedronid = dao_polyhedron->createPolyhedronID(facet->getPolyhedron());
    }
    if (polyhedronid > 0) {
        int fid = nextFID(polyhedronid);
        std::string sql("INSERT INTO Facets (PolyhedronID, FID) "
                "VALUES (?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, polyhedronid);
            stmt->bindInteger(2, fid);
            if (stmt->execute()) {
                result = fid;
                facet->setID(result);
            }
        }
    }
    return result;
}

int FacetDAO::insert(FacetSPtr facet) {
    int result = -1;
    if (!facet->getPolyhedron()) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int fid = createFID(facet);
    if (fid > 0) {
        VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
        std::list<VertexSPtr>::iterator it_v = facet->vertices().begin();
        while (it_v != facet->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->getID() > 0) {
                dao_vertex->update(vertex);
            } else {
                dao_vertex->insert(vertex);
            }
        }
        EdgeDAOSPtr dao_edge = DAOFactory::getEdgeDAO();
        std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
        while (it_e != facet->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (edge->getID() > 0) {
                dao_edge->update(edge);
            } else {
                dao_edge->insert(edge);
            }
        }
        TriangleDAOSPtr dao_triangle = DAOFactory::getTriangleDAO();
        std::list<TriangleSPtr>::iterator it_t = facet->triangles().begin();
        while (it_t != facet->triangles().end()) {
            TriangleSPtr triangle = *it_t++;
            if (triangle->getID() > 0) {
                dao_triangle->update(triangle);
            } else {
                dao_triangle->insert(triangle);
            }
        }
        Plane3SPtr plane = facet->getPlane();
        if (plane) {
            PlaneDAOSPtr dao_plane = DAOFactory::getPlaneDAO();
            int plane_id = dao_plane->insert(plane);
            if (plane_id > 0) {
                std::string sql("UPDATE Facets SET PlaneID=? "
                        "WHERE PolyhedronID=? AND FID=?;");
                SQLiteStmtSPtr stmt = db->prepare(sql);
                if (stmt) {
                    stmt->bindInteger(1, plane_id);
                    stmt->bindInteger(2, facet->getPolyhedron()->getID());
                    stmt->bindInteger(3, facet->getID());
                    stmt->execute();
                }
            }
        }
        result = fid;
    }
    return result;
}

bool FacetDAO::del(FacetSPtr facet) {
    bool result = false;
    if (!facet->getPolyhedron()) {
        return false;
    }
    int polyhedronid = facet->getPolyhedron()->getID();
    int fid = facet->getID();
    if (polyhedronid < 0 || fid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("DELETE FROM Facets WHERE PolyhedronID=? AND FID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, fid);
        if (stmt->execute() > 0) {
            sql = "UPDATE Edges SET FID_L=? "
                    "WHERE PolyhedronID=? AND FID_L=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, 0);
                stmt->bindInteger(2, polyhedronid);
                stmt->bindInteger(3, fid);
                stmt->execute();
            }
            sql = "UPDATE Edges SET FID_R=? "
                    "WHERE PolyhedronID=? AND FID_R=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, 0);
                stmt->bindInteger(2, polyhedronid);
                stmt->bindInteger(3, fid);
                stmt->execute();
            }
            sql = "DELETE FROM Triangles WHERE PolyhedronID=? AND "
                    "(FID_L=? OR FID_R=?);";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedronid);
                stmt->bindInteger(2, fid);
                stmt->bindInteger(3, fid);
                stmt->execute();
            }
            result = true;
            facet->setID(-1);
        }
    }
    return result;
}

FacetSPtr FacetDAO::find(int polyhedronid, int fid) {
    FacetSPtr result = FacetSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
    std::string sql("SELECT FID, PlaneID FROM Facets "
            "WHERE PolyhedronID=? AND FID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        stmt->bindInteger(2, fid);
        if (stmt->execute() > 0) {
            int plane_id = stmt->getInteger(1);
            result = Facet::create();
            result->setID(fid);
            if (plane_id > 0) {
                PlaneDAOSPtr dao_plane = DAOFactory::getPlaneDAO();
                Plane3SPtr plane = dao_plane->find(plane_id);
                result->setPlane(plane);
            }
            std::map<int, VertexSPtr> vertices;
            sql = "SELECT EID, VID_SRC, VID_DST, FID_L, FID_R FROM Edges "
                    "WHERE PolyhedronID=? AND (FID_L=? OR FID_R=?);";
            SQLiteStmtSPtr stmt_e = db->prepare(sql);
            if (stmt_e) {
                stmt_e->bindInteger(1, polyhedronid);
                stmt_e->bindInteger(2, polyhedronid);
                stmt_e->bindInteger(3, polyhedronid);
                stmt_e->execute();
                while (stmt_e->fetchRow()) {
                    int eid = stmt_e->getInteger(0);
                    int vid_src = stmt_e->getInteger(1);
                    int vid_dst = stmt_e->getInteger(2);
                    int fid_l = stmt_e->getInteger(3);
                    int fid_r = stmt_e->getInteger(4);
                    if (!vertices[vid_src]) {
                        VertexSPtr vertex = dao_vertex->find(polyhedronid, vid_src);
                        vertices[vid_src] = vertex;
                        result->addVertex(vertex);
                    }
                    if (!vertices[vid_dst]) {
                        VertexSPtr vertex = dao_vertex->find(polyhedronid, vid_dst);
                        vertices[vid_dst] = vertex;
                        result->addVertex(vertex);
                    }
                    EdgeSPtr edge = Edge::create(vertices[vid_src], vertices[vid_dst]);
                    edge->setID(eid);
                    if (fid_l == fid) {
                        edge->setFacetL(result);
                    }
                    if (fid_r == fid) {
                        edge->setFacetR(result);
                    }
                    result->addEdge(edge);
                }
                stmt_e->close();
            }
            sql = "SELECT TID, VID_1, VID_2, VID_3 FROM Triangles "
                    "WHERE PolyhedronID=? AND FID=?;";
            SQLiteStmtSPtr stmt_t = db->prepare(sql);
            if (stmt_t) {
                stmt_t->bindInteger(1, polyhedronid);
                stmt_t->bindInteger(2, fid);
                stmt_t->execute();
                while (stmt_t->fetchRow()) {
                    int tid = stmt_t->getInteger(0);
                    int vid1 = stmt_t->getInteger(1);
                    int vid2 = stmt_t->getInteger(2);
                    int vid3 = stmt_t->getInteger(3);
                    VertexSPtr verts[3];
                    verts[0] = vertices[vid1];
                    verts[1] = vertices[vid2];
                    verts[2] = vertices[vid3];
                    TriangleSPtr triangle = Triangle::create(result, verts);
                    triangle->setID(tid);
                }
            }
            result->makeFirstConvex();
        }
    }
    return result;
}

bool FacetDAO::update(FacetSPtr facet) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
