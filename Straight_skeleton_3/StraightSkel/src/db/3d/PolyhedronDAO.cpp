#include "db/3d/PolyhedronDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/VertexDAO.h"
#include "db/3d/EdgeDAO.h"
#include "db/3d/FacetDAO.h"
#include "db/3d/TriangleDAO.h"
#include "db/3d/PlaneDAO.h"
#include <map>

namespace db { namespace _3d {

PolyhedronDAO::PolyhedronDAO() {
    // intentionally does nothing
}

PolyhedronDAO::~PolyhedronDAO() {
    // intentionally does nothing
}

std::string PolyhedronDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Polyhedrons (\n"
            "  PolyhedronID INTEGER PRIMARY KEY,\n"
            "  description TEXT,\n"
            "  created INTEGER\n"
            ");");
    return schema;
}

int PolyhedronDAO::nextPolyhedronID() {
    int polyhedronid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(PolyhedronID) FROM Polyhedrons;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        polyhedronid = 1;
        if (stmt->execute()) {
            polyhedronid = stmt->getInteger(0) + 1;
        }
    }
    return polyhedronid;
}

int PolyhedronDAO::createPolyhedronID(PolyhedronSPtr polyhedron) {
    int result = -1;
    int polyhedronid = nextPolyhedronID();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("INSERT INTO Polyhedrons (PolyhedronID, description, created) "
            "VALUES (?, ?, strftime('%s','now'));");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        if (polyhedronid > 0) {
            stmt->bindInteger(1, polyhedronid);
            stmt->bindString(2, polyhedron->getDescription());
            if (stmt->execute()) {
                result = polyhedronid;
                polyhedron->setID(result);
            }
        }
    }
    return result;
}

int PolyhedronDAO::insert(PolyhedronSPtr polyhedron) {
    int result = -1;
    WriteLock l(polyhedron->mutex());
    polyhedron->resetAllIDs();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    int polyhedronid = createPolyhedronID(polyhedron);
    if (polyhedronid > 0) {
        VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
        std::list<VertexSPtr>::iterator it_v = polyhedron->vertices().begin();
        while (it_v != polyhedron->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->getID() > 0) {
                dao_vertex->update(vertex);
            } else {
                dao_vertex->insert(vertex);
            }
        }
        EdgeDAOSPtr dao_edge = DAOFactory::getEdgeDAO();
        std::list<EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
        while (it_e != polyhedron->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (edge->getID() > 0) {
                dao_edge->update(edge);
            } else {
                dao_edge->insert(edge);
            }
        }
        FacetDAOSPtr dao_facet = DAOFactory::getFacetDAO();
        std::list<FacetSPtr>::iterator it_f = polyhedron->facets().begin();
        while (it_f != polyhedron->facets().end()) {
            FacetSPtr facet = *it_f++;
            if (facet->getID() > 0) {
                dao_facet->update(facet);
            } else {
                dao_facet->insert(facet);
            }
        }
        result = polyhedronid;
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

bool PolyhedronDAO::del(PolyhedronSPtr polyhedron) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    std::string sql("DELETE FROM Polyhedrons WHERE PolyhedronID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedron->getID());
        if (stmt->execute() > 0) {
            sql = "DELETE FROM Edges WHERE PolyhedronID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedron->getID());
                stmt->execute();
            }
            sql = "DELETE FROM Triangles WHERE PolyhedronID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedron->getID());
                stmt->execute();
            }
            sql = "DELETE FROM Facets WHERE PolyhedronID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedron->getID());
                stmt->execute();
            }
            sql = "DELETE FROM Vertices WHERE PolyhedronID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyhedron->getID());
                stmt->execute();
            }
            result = true;
            polyhedron->setID(-1);
        }
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

PolyhedronSPtr PolyhedronDAO::find(int polyhedronid) {
    PolyhedronSPtr result = PolyhedronSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    PlaneDAOSPtr plane_dao = DAOFactory::getPlaneDAO();
    VertexDAOSPtr vertex_dao = DAOFactory::getVertexDAO();
    std::string sql("SELECT PolyhedronID, description FROM Polyhedrons WHERE PolyhedronID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyhedronid);
        if (stmt->execute() > 0) {
            result = Polyhedron::create();
            result->setID(polyhedronid);
            std::string description = stmt->getString(1);
            result->setDescription(description);
            std::map<int, VertexSPtr> vertices;
            vertices.clear();
            sql = "SELECT VID FROM Vertices WHERE PolyhedronID=? ORDER BY VID ASC;";
            SQLiteStmtSPtr stmt_v = db->prepare(sql);
            if (stmt_v) {
                stmt_v->bindInteger(1, polyhedronid);
                stmt_v->execute();
                while (stmt_v->fetchRow()) {
                    int vid = stmt_v->getInteger(0);
                    VertexSPtr vertex = vertex_dao->find(polyhedronid, vid);
                    vertices[vid] = vertex;
                    result->addVertex(vertex);
                }
                stmt_v->close();
            }
            std::map<int, FacetSPtr> facets;
            sql = "SELECT FID, PlaneID FROM Facets WHERE PolyhedronID=?;";
            SQLiteStmtSPtr stmt_p = db->prepare(sql);
            if (stmt_p) {
                stmt_p->bindInteger(1, polyhedronid);
                stmt_p->execute();
                while (stmt_p->fetchRow()) {
                    int fid = stmt_p->getInteger(0);
                    int plane_id = stmt_p->getInteger(1);
                    FacetSPtr facet = Facet::create();
                    facet->setID(fid);
                    if (plane_id > 0) {
                        Plane3SPtr plane = plane_dao->find(plane_id);
                        facet->setPlane(plane);
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
                            int vids[3];
                            for (unsigned int i = 0; i < 3; i++) {
                                vids[i] = stmt_t->getInteger(i+1);
                            }
                            VertexSPtr verts[3];
                            for (unsigned int i = 0; i < 3; i++) {
                                verts[i] = vertices[vids[i]];
                            }
                            TriangleSPtr triangle = Triangle::create(facet, verts);
                            triangle->setID(tid);
                        }
                        stmt_t->close();
                    }
                    facets[fid] = facet;
                    result->addFacet(facet);
                }
                stmt_p->close();
            }
            sql = "SELECT EID, VID_SRC, VID_DST, FID_L, FID_R "
                    "FROM Edges WHERE PolyhedronID=? ORDER BY EID ASC;";
            SQLiteStmtSPtr stmt_e = db->prepare(sql);
            if (stmt_e) {
                stmt_e->bindInteger(1, polyhedronid);
                stmt_e->execute();
                while (stmt_e->fetchRow()) {
                    int eid = stmt_e->getInteger(0);
                    int vid_src = stmt_e->getInteger(1);
                    int vid_dst = stmt_e->getInteger(2);
                    int fid_l = stmt_e->getInteger(3);
                    int fid_r = stmt_e->getInteger(4);
                    EdgeSPtr edge = Edge::create(vertices[vid_src], vertices[vid_dst]);
                    if (fid_l > 0) {
                        edge->setFacetL(facets[fid_l]);
                        edge->getFacetL()->addEdge(edge);
                    }
                    if (fid_r > 0) {
                        edge->setFacetR(facets[fid_r]);
                        edge->getFacetR()->addEdge(edge);
                    }
                    result->addEdge(edge);
                    edge->setID(eid);
                }
                stmt_e->close();
            }
            std::list<FacetSPtr>::iterator it_f = result->facets().begin();
            while (it_f != result->facets().end()) {
                FacetSPtr facet = *it_f++;
                facet->makeFirstConvex();
            }
        }
    }
    return result;
}

bool PolyhedronDAO::update(PolyhedronSPtr polyhedron) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
