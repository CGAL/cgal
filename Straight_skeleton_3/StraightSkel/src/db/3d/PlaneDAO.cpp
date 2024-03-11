/**
 * @file   db/3d/PlaneDAO.cpp
 * @author Gernot Walzl
 * @date   2013-06-07
 */

#include "db/3d/PlaneDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"

namespace db { namespace _3d {

PlaneDAO::PlaneDAO() {
    // intentionally does nothing
}

PlaneDAO::~PlaneDAO() {
    // intentionally does nothing
}

std::string PlaneDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Planes (\n"
            "  PlaneID INTEGER PRIMARY KEY,\n"
            "  a REAL,\n"
            "  b REAL,\n"
            "  c REAL,\n"
            "  d REAL\n"
            ");");
    return schema;
}

int PlaneDAO::nextPlaneID() {
    int plane_id = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(PlaneID) FROM Planes;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        plane_id = 1;
        if (stmt->execute()) {
            plane_id = stmt->getInteger(0) + 1;
        }
    }
    return plane_id;
}

int PlaneDAO::insert(Plane3SPtr plane) {
    int result = -1;
    if (!plane) {
        return -1;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int plane_id = nextPlaneID();
    if (plane_id > 0) {
        std::string sql("INSERT INTO Planes (PlaneID, a, b, c, d) "
                "VALUES (?, ?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            double d = 0.0;
#ifdef USE_CGAL
            a = plane->a();
            b = plane->b();
            c = plane->c();
            d = plane->d();
#else
            a = plane->getA();
            b = plane->getB();
            c = plane->getC();
            d = plane->getD();
#endif
            stmt->bindInteger(1, plane_id);
            stmt->bindDouble(2, a);
            stmt->bindDouble(3, b);
            stmt->bindDouble(4, c);
            stmt->bindDouble(5, d);
            if (stmt->execute() > 0) {
                result = plane_id;
            }
        }
    }
    return result;
}

bool PlaneDAO::del(Plane3SPtr plane) {
    bool result = false;
    int plane_id = 0;
    // TODO
    if (plane_id > 0) {
        SQLiteDatabaseSPtr db = DAOFactory::getDB();
        std::string sql("DELETE FROM Planes WHERE PlaneID=?;");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, plane_id);
            if (stmt->execute() > 0) {
                result = true;
            }
        }
    }
    return result;
}

Plane3SPtr PlaneDAO::find(int plane_id) {
    Plane3SPtr result = Plane3SPtr();
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT a, b, c, d FROM Planes WHERE PlaneID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, plane_id);
        if (stmt->execute() > 0) {
            double a = stmt->getDouble(0);
            double b = stmt->getDouble(1);
            double c = stmt->getDouble(2);
            double d = stmt->getDouble(3);
            result = KernelFactory::createPlane3(a,b,c,d);
        }
    }
    return result;
}

bool PlaneDAO::update(Plane3SPtr point) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
