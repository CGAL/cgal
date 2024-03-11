/**
 * @file   db/3d/DAOFactory.cpp
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#include "db/3d/DAOFactory.h"

#include "db/SQLiteDatabase.h"
#include "db/SQLiteStmt.h"
#include "db/3d/PointDAO.h"
#include "db/3d/PlaneDAO.h"
#include "db/3d/VertexDAO.h"
#include "db/3d/EdgeDAO.h"
#include "db/3d/TriangleDAO.h"
#include "db/3d/FacetDAO.h"
#include "db/3d/PolyhedronDAO.h"
#include "db/3d/NodeDAO.h"
#include "db/3d/ArcDAO.h"
#include "db/3d/SheetDAO.h"
#include "db/3d/EventDAO.h"
#include "db/3d/StraightSkeletonDAO.h"
#include <cstdlib>
#include <fstream>

namespace db { namespace _3d {

SQLiteDatabaseSPtr DAOFactory::db_;

PointDAOSPtr DAOFactory::point_dao_;
PlaneDAOSPtr DAOFactory::plane_dao_;

VertexDAOSPtr DAOFactory::vertex_dao_;
EdgeDAOSPtr DAOFactory::edge_dao_;
TriangleDAOSPtr DAOFactory::triangle_dao_;
FacetDAOSPtr DAOFactory::facet_dao_;
PolyhedronDAOSPtr DAOFactory::polyhedron_dao_;

NodeDAOSPtr DAOFactory::node_dao_;
ArcDAOSPtr DAOFactory::arc_dao_;
SheetDAOSPtr DAOFactory::sheet_dao_;
EventDAOSPtr DAOFactory::event_dao_;
StraightSkeletonDAOSPtr DAOFactory::straightskeleton_dao_;


DAOFactory::DAOFactory() {
    // intentionally does nothing
}

DAOFactory::~DAOFactory() {
    straightskeleton_dao_.reset();
    event_dao_.reset();
    sheet_dao_.reset();
    arc_dao_.reset();
    node_dao_.reset();

    polyhedron_dao_.reset();
    facet_dao_.reset();
    triangle_dao_.reset();
    edge_dao_.reset();
    vertex_dao_.reset();

    plane_dao_.reset();
    point_dao_.reset();

    db_.reset();
}

std::string DAOFactory::findDefaultFilename() {
    std::string name("StraightSkel");
    std::string filename("skeldata3d.db3");
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
    std::string sql[13];
    sql[0] = getPointDAO()->getTableSchema();
    sql[1] = getPlaneDAO()->getTableSchema();
    sql[2] = getVertexDAO()->getTableSchema();
    sql[3] = getEdgeDAO()->getTableSchema();
    sql[4] = getTriangleDAO()->getTableSchema();
    sql[5] = getFacetDAO()->getTableSchema();
    sql[6] = getPolyhedronDAO()->getTableSchema();
    sql[7] = getNodeDAO()->getTableSchema();
    sql[8] = getArcDAO()->getTableSchema();
    sql[9] = getSheetDAO()->getTableSchema();
    sql[10] = getSheetDAO()->getTable2Schema();
    sql[11] = getEventDAO()->getTableSchema();
    sql[12] = getStraightSkeletonDAO()->getTableSchema();
    db_->beginTransaction();
    for (unsigned int i = 0; i < 13; i++) {
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

PlaneDAOSPtr DAOFactory::getPlaneDAO() {
    if (!plane_dao_) {
        plane_dao_ = PlaneDAOSPtr(new PlaneDAO());
    }
    return plane_dao_;
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

TriangleDAOSPtr DAOFactory::getTriangleDAO() {
    if (!triangle_dao_) {
       triangle_dao_ = TriangleDAOSPtr(new TriangleDAO());
    }
    return triangle_dao_;
}

FacetDAOSPtr DAOFactory::getFacetDAO() {
    if (!facet_dao_) {
        facet_dao_ = FacetDAOSPtr(new FacetDAO());
    }
    return facet_dao_;
}

PolyhedronDAOSPtr DAOFactory::getPolyhedronDAO() {
  if (!polyhedron_dao_) {
        polyhedron_dao_ = PolyhedronDAOSPtr(new PolyhedronDAO());
    }
    return polyhedron_dao_;
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

SheetDAOSPtr DAOFactory::getSheetDAO() {
    if (!sheet_dao_) {
        sheet_dao_ = SheetDAOSPtr(new SheetDAO());
    }
    return sheet_dao_;
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
