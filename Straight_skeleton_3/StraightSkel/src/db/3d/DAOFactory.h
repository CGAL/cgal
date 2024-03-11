/**
 * @file   db/3d/DAOFactory.h
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#ifndef DB_3D_DAOFACTORY_H
#define DB_3D_DAOFACTORY_H


#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include <string>

namespace db { namespace _3d {

class DAOFactory {
public:
    virtual ~DAOFactory();

    static std::string findDefaultFilename();
    static bool createTables();
    static SQLiteDatabaseSPtr getDB();

    static PointDAOSPtr getPointDAO();
    static PlaneDAOSPtr getPlaneDAO();

    static VertexDAOSPtr getVertexDAO();
    static EdgeDAOSPtr getEdgeDAO();
    static TriangleDAOSPtr getTriangleDAO();
    static FacetDAOSPtr getFacetDAO();
    static PolyhedronDAOSPtr getPolyhedronDAO();

    static NodeDAOSPtr getNodeDAO();
    static ArcDAOSPtr getArcDAO();
    static SheetDAOSPtr getSheetDAO();
    static EventDAOSPtr getEventDAO();
    static StraightSkeletonDAOSPtr getStraightSkeletonDAO();

protected:
    DAOFactory();

    static SQLiteDatabaseSPtr db_;

    static PointDAOSPtr point_dao_;
    static PlaneDAOSPtr plane_dao_;

    static VertexDAOSPtr vertex_dao_;
    static EdgeDAOSPtr edge_dao_;
    static TriangleDAOSPtr triangle_dao_;
    static FacetDAOSPtr facet_dao_;
    static PolyhedronDAOSPtr polyhedron_dao_;

    static NodeDAOSPtr node_dao_;
    static ArcDAOSPtr arc_dao_;
    static SheetDAOSPtr sheet_dao_;
    static EventDAOSPtr event_dao_;
    static StraightSkeletonDAOSPtr straightskeleton_dao_;
};

} }

#endif /* DB_3D_DAOFACTORY_H */
