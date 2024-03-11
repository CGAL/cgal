/**
 * @file   db/3d/PointDAO.h
 * @author Gernot Walzl
 * @date   2013-05-24
 */

#ifndef DB_3D_POINTDAO_H
#define DB_3D_POINTDAO_H

#include <string>
#include <map>
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"

namespace db { namespace _3d {

using data::_3d::KernelFactory;
using data::_3d::Point3SPtr;
using data::_3d::Vector3SPtr;

class PointDAO {
friend class DAOFactory;
public:
    virtual ~PointDAO();
    std::string getTableSchema() const;
    int insert(Point3SPtr point);
    bool del(Point3SPtr point);
    Point3SPtr find(int point_id);
    bool update(Point3SPtr point);
protected:
    PointDAO();
    int nextPointID();
    static std::map<Point3SPtr, int> point_ids_;
    static std::map<int, Point3SPtr> points_;
};

} }

#endif /* DB_3D_POINTDAO_H */
