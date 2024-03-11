/**
 * @file   db/3d/PlaneDAO.h
 * @author Gernot Walzl
 * @date   2013-06-07
 */

#ifndef DB_3D_PLANEDAO_H
#define DB_3D_PLANEDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::KernelFactory;
using data::_3d::Plane3SPtr;

class PlaneDAO {
friend class DAOFactory;
public:
    virtual ~PlaneDAO();
    std::string getTableSchema() const;
    int insert(Plane3SPtr plane);
    bool del(Plane3SPtr plane);
    Plane3SPtr find(int plane_id);
    bool update(Plane3SPtr plane);
protected:
    PlaneDAO();
    int nextPlaneID();
};

} }

#endif /* DB_3D_PLANEDAO_H */

