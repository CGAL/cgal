/**
 * @file   db/2d/PolygonDAO.cpp
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#include "db/2d/PolygonDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/2d/VertexDAO.h"
#include "db/2d/EdgeDAO.h"
#include <list>
#include <map>

namespace db { namespace _2d {

PolygonDAO::PolygonDAO() {
    // intentionally does nothing
}

PolygonDAO::~PolygonDAO() {
    // intentionally does nothing
}

std::string PolygonDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Polygons (\n"
            "  PolyID INTEGER PRIMARY KEY,\n"
            "  description TEXT,\n"
            "  created INTEGER\n"
            ");");
    return schema;
}

int PolygonDAO::nextPolyID() {
    int polyid = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(PolyID) FROM Polygons;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        polyid = 1;
        if (stmt->execute()) {
            polyid = stmt->getInteger(0) + 1;
        }
    }
    return polyid;
}

int PolygonDAO::createPolyID(PolygonSPtr polygon) {
    int result = -1;
    int polyid = nextPolyID();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("INSERT INTO Polygons (PolyID, description, created) "
            "VALUES (?, ?, strftime('%s','now'));");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        if (polyid > 0) {
            stmt->bindInteger(1, polyid);
            stmt->bindString(2, polygon->getDescription());
            if (stmt->execute()) {
                result = polyid;
                polygon->setID(result);
            }
        }
    }
    return result;
}

int PolygonDAO::insert(PolygonSPtr polygon) {
    int result = -1;
    WriteLock(polygon->mutex());
    polygon->resetAllIDs();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    int polyid = createPolyID(polygon);
    if (polyid > 0) {
        VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
        std::list<VertexSPtr>::iterator it_v = polygon->vertices().begin();
        while (it_v != polygon->vertices().end()) {
            VertexSPtr vertex = *it_v++;
            if (vertex->getID() > 0) {
                dao_vertex->update(vertex);
            } else {
                dao_vertex->insert(vertex);
            }
        }
        EdgeDAOSPtr dao_edge = DAOFactory::getEdgeDAO();
        std::list<EdgeSPtr>::iterator it_e = polygon->edges().begin();
        while (it_e != polygon->edges().end()) {
            EdgeSPtr edge = *it_e++;
            if (edge->getID() > 0) {
                dao_edge->update(edge);
            } else {
                dao_edge->insert(edge);
            }
        }
        result = polyid;
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

bool PolygonDAO::del(PolygonSPtr polygon) {
    bool result = false;
    int polyid = polygon->getID();
    if (polyid < 0) {
        return false;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    bool trans_started = db->beginTransaction();
    std::string sql("DELETE FROM Polygons WHERE PolyID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        if (stmt->execute() > 0) {
            sql = "DELETE FROM Edges WHERE PolyID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyid);
                stmt->execute();
            }
            sql = "DELETE FROM Vertices WHERE PolyID=?;";
            stmt = db->prepare(sql);
            if (stmt) {
                stmt->bindInteger(1, polyid);
                stmt->execute();
            }
            result = true;
            polygon->setID(-1);
        }
    }
    if (trans_started) {
        db->endTransaction();
    }
    return result;
}

PolygonSPtr PolygonDAO::find(int polyid) {
    PolygonSPtr result = PolygonSPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    VertexDAOSPtr dao_vertex = DAOFactory::getVertexDAO();
    std::string sql("SELECT PolyID, description FROM Polygons WHERE PolyID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, polyid);
        if (stmt->execute() > 0) {
            result = Polygon::create();
            result->setID(polyid);
            std::string description = stmt->getString(1);
            result->setDescription(description);
            std::map<int, VertexSPtr> vertices;
            sql = "SELECT VID FROM Vertices WHERE PolyID=? ORDER BY VID ASC;";
            SQLiteStmtSPtr stmt_v = db->prepare(sql);
            if (stmt_v) {
                stmt_v->bindInteger(1, polyid);
                stmt_v->execute();
                while (stmt_v->fetchRow()) {
                    int vid = stmt_v->getInteger(0);
                    VertexSPtr vertex = dao_vertex->find(polyid, vid);
                    vertices[vid] = vertex;
                    result->addVertex(vertex);
                }
                stmt_v->close();
            }
            sql = "SELECT EID, VID_SRC, VID_DST FROM Edges WHERE PolyID=? ORDER BY EID ASC;";
            SQLiteStmtSPtr stmt_e = db->prepare(sql);
            if (stmt_e) {
                stmt_e->bindInteger(1, polyid);
                stmt_e->execute();
                while (stmt_e->fetchRow()) {
                    int eid = stmt_e->getInteger(0);
                    int vid_src = stmt_e->getInteger(1);
                    int vid_dst = stmt_e->getInteger(2);
                    EdgeSPtr edge = Edge::create(vertices[vid_src], vertices[vid_dst]);
                    edge->setID(eid);
                    result->addEdge(edge);
                    sql = "SELECT speed FROM SkelEdgeData WHERE PolyID=? AND EID=?;";
                    SQLiteStmtSPtr stmt_d = db->prepare(sql);
                    if (stmt_d) {
                        stmt_d->bindInteger(1, polyid);
                        stmt_d->bindInteger(2, eid);
                        if (stmt_d->execute() > 0) {
                            double speed = stmt_d->getDouble(0);
                            if (speed != 0.0 && speed != 1.0) {
                                SkelEdgeDataSPtr data = SkelEdgeData::create(edge);
                                data->setSpeed(speed);
                            }
                        }
                        stmt_d->close();
                    }
                }
                stmt_e->close();
            }
        }
    }
    return result;
}

bool PolygonDAO::update(PolygonSPtr polygon) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
