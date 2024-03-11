/**
 * @file   db/2d/DAOFactory.cpp
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#include "db/2d/DAOFactory.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/2d/PointDAO.h"
#include "db/2d/VertexDAO.h"
#include "db/2d/EdgeDAO.h"
#include "db/2d/PolygonDAO.h"
#include "db/2d/NodeDAO.h"
#include "db/2d/ArcDAO.h"
#include "db/2d/EventDAO.h"
#include "db/2d/StraightSkeletonDAO.h"
#include <cstdlib>
#include <fstream>

namespace db { namespace _2d {

SQLiteDatabaseSPtr DAOFactory::db_;

PointDAOSPtr DAOFactory::point_dao_;

VertexDAOSPtr DAOFactory::vertex_dao_;
EdgeDAOSPtr DAOFactory::edge_dao_;
PolygonDAOSPtr DAOFactory::polygon_dao_;

NodeDAOSPtr DAOFactory::node_dao_;
ArcDAOSPtr DAOFactory::arc_dao_;
EventDAOSPtr DAOFactory::event_dao_;
StraightSkeletonDAOSPtr DAOFactory::straightskeleton_dao_;


DAOFactory::DAOFactory() {
    // intentionally does nothing
}

DAOFactory::~DAOFactory() {
    straightskeleton_dao_.reset();
    event_dao_.reset();
    arc_dao_.reset();
    node_dao_.reset();

    polygon_dao_.reset();
    edge_dao_.reset();
    vertex_dao_.reset();

    point_dao_.reset();

    db_.reset();
}

std::string DAOFactory::findDefaultFilename() {
    std::string name("StraightSkel");
    std::string filename("skeldata2d.db3");
    std::string result = filename;
    std::string home(getenv("HOME"));
//    boost::filesystem::path p(home+"/."+name);
//    if (boost::filesystem::is_directory(p)) {
//        result = home+"/."+name+"/"+filename;
//    }
    std::string filenames[2];
    filenames[0] = filename;
    filenames[1] = home+"/."+name+"/"+filename;
    for (unsigned int i = 0; i < 1; i++) {
        std::ifstream input(filenames[i].c_str());
        if (input.is_open()) {
            input.close();
            result = filenames[i];
            break;
        }
    }
    return result;
}

bool DAOFactory::createTables() {
    bool result = true;
    std::string sql[9];
    sql[0] = getPointDAO()->getTableSchema();
    sql[1] = getVertexDAO()->getTableSchema();
    sql[2] = getEdgeDAO()->getTableSchema();
    sql[3] = getEdgeDAO()->getTable2Schema();
    sql[4] = getPolygonDAO()->getTableSchema();
    sql[5] = getNodeDAO()->getTableSchema();
    sql[6] = getArcDAO()->getTableSchema();
    sql[7] = getEventDAO()->getTableSchema();
    sql[8] = getStraightSkeletonDAO()->getTableSchema();
    db_->beginTransaction();
    for (unsigned int i = 0; i < 9; i++) {
        SQLiteStmtSPtr stmt = db_->prepare(sql[i]);
        if (!stmt->execute()) {
            result = false;
        }
    }
    db_->endTransaction();
    return result;
}

SQLiteDatabaseSPtr DAOFactory::getDB() {
    if (!db_) {
        db_ = SQLiteDatabaseSPtr(new SQLiteDatabase());
        std::string filename = findDefaultFilename();
        std::ifstream input(filename.c_str());
        if (input.is_open()) {
            input.close();
            db_->open(filename);
        } else {
            db_->open(filename);
            createTables();
        }
    }
    return db_;
}

PointDAOSPtr DAOFactory::getPointDAO() {
    if (!point_dao_) {
        point_dao_ = PointDAOSPtr(new PointDAO());
    }
    return point_dao_;
}

VertexDAOSPtr DAOFactory::getVertexDAO() {
    if (!vertex_dao_) {
        vertex_dao_ = VertexDAOSPtr(new VertexDAO());
    }
    return vertex_dao_;
}

EdgeDAOSPtr DAOFactory::getEdgeDAO() {
    if (!edge_dao_) {
        edge_dao_ = EdgeDAOSPtr(new EdgeDAO());
    }
    return edge_dao_;
}

PolygonDAOSPtr DAOFactory::getPolygonDAO() {
    if (!polygon_dao_) {
        polygon_dao_ = PolygonDAOSPtr(new PolygonDAO());
    }
    return polygon_dao_;
}

NodeDAOSPtr DAOFactory::getNodeDAO() {
    if (!node_dao_) {
        node_dao_ = NodeDAOSPtr(new NodeDAO());
    }
    return node_dao_;
}

ArcDAOSPtr DAOFactory::getArcDAO() {
    if (!arc_dao_) {
        arc_dao_ = ArcDAOSPtr(new ArcDAO());
    }
    return arc_dao_;
}

EventDAOSPtr DAOFactory::getEventDAO() {
    if (!event_dao_) {
        event_dao_ = EventDAOSPtr(new EventDAO());
    }
    return event_dao_;
}

StraightSkeletonDAOSPtr DAOFactory::getStraightSkeletonDAO() {
    if (!straightskeleton_dao_) {
        straightskeleton_dao_ = StraightSkeletonDAOSPtr(new StraightSkeletonDAO());
    }
    return straightskeleton_dao_;
}

} }
