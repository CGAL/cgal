/**
 * @file   db/3d/PointDAO.cpp
 * @author Gernot Walzl
 * @date   2013-05-24
 */

#include "db/3d/PointDAO.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"

namespace db { namespace _3d {

std::map<Point3SPtr, int> PointDAO::point_ids_;
std::map<int, Point3SPtr> PointDAO::points_;

PointDAO::PointDAO() {
    // intentionally does nothing
}

PointDAO::~PointDAO() {
    point_ids_.clear();
    points_.clear();
}

std::string PointDAO::getTableSchema() const {
    std::string schema("CREATE TABLE Points (\n"
            "  PointID INTEGER PRIMARY KEY,\n"
            "  x REAL,\n"
            "  y REAL,\n"
            "  z REAL\n"
            ");");
    return schema;
}

int PointDAO::nextPointID() {
    int point_id = -1;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT MAX(PointID) FROM Points;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        point_id = 1;
        if (stmt->execute()) {
            point_id = stmt->getInteger(0) + 1;
        }
    }
    return point_id;
}

int PointDAO::insert(Point3SPtr point) {
    int result = -1;
    if (!point) {
        return -1;
    }
    std::map<Point3SPtr, int>::const_iterator it_p = point_ids_.find(point);
    if (it_p != point_ids_.end()) {
        result = it_p->second;
        return result;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    int point_id = nextPointID();
    if (point_id > 0) {
        std::string sql("INSERT INTO Points (PointID, x, y, z) VALUES (?, ?, ?, ?);");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            Vector3SPtr vec = KernelFactory::createVector3(point);
            stmt->bindInteger(1, point_id);
            stmt->bindDouble(2, (*vec)[0]);
            stmt->bindDouble(3, (*vec)[1]);
            stmt->bindDouble(4, (*vec)[2]);
            if (stmt->execute() > 0) {
                point_ids_[point] = point_id;
                result = point_id;
            }
        }
    }
    return result;
}

bool PointDAO::del(Point3SPtr point) {
    bool result = false;
    int point_id = 0;
    std::map<Point3SPtr, int>::iterator it_p = point_ids_.find(point);
    if (it_p != point_ids_.end()) {
        point_id = it_p->second;
    }
    if (point_id > 0) {
        SQLiteDatabaseSPtr db = DAOFactory::getDB();
        std::string sql("DELETE FROM Points WHERE PointID=?;");
        SQLiteStmtSPtr stmt = db->prepare(sql);
        if (stmt) {
            stmt->bindInteger(1, point_id);
            if (stmt->execute() > 0) {
                points_.erase(points_.find(point_id));
                point_ids_.erase(it_p);
                result = true;
            }
        }
    }
    return result;
}

Point3SPtr PointDAO::find(int point_id) {
    Point3SPtr result = Point3SPtr();
    std::map<int, Point3SPtr>::const_iterator it_p = points_.find(point_id);
    if (it_p != points_.end()) {
        result = it_p->second;
        return result;
    }
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    std::string sql("SELECT x, y, z FROM Points WHERE PointID=?;");
    SQLiteStmtSPtr stmt = db->prepare(sql);
    if (stmt) {
        stmt->bindInteger(1, point_id);
        if (stmt->execute() > 0) {
            double x = stmt->getDouble(0);
            double y = stmt->getDouble(1);
            double z = stmt->getDouble(2);
            result = KernelFactory::createPoint3(x,y,z);
            points_[point_id] = result;
        }
    }
    return result;
}

bool PointDAO::update(Point3SPtr point) {
    bool result = false;
    SQLiteDatabaseSPtr db = DAOFactory::getDB();
    // TODO
    return result;
}

} }
