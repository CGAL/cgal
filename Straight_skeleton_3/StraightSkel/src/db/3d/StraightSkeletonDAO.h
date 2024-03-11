/**
 * @file   db/3d/StraightSkeletonDAO.h
 * @author Gernot Walzl
 * @date   2013-06-03
 */

#ifndef DB_3D_STRAIGHTSKELETONDAO_H
#define DB_3D_STRAIGHTSKELETONDAO_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"

namespace db { namespace _3d {

using data::_3d::skel::StraightSkeleton;
using data::_3d::skel::StraightSkeletonSPtr;

class StraightSkeletonDAO {
friend class DAOFactory;
public:
    virtual ~StraightSkeletonDAO();
    std::string getTableSchema() const;
    int createSkelID(StraightSkeletonSPtr skel);
    int insert(StraightSkeletonSPtr skel);
    bool del(StraightSkeletonSPtr skel);
    StraightSkeletonSPtr find(int skelid);
    int findPolyhedronID(int skelid);
    bool update(StraightSkeletonSPtr skel);
private:
    StraightSkeletonDAO();
    int nextSkelID();
};

} }

#endif /* DB_3D_STRAIGHTSKELETONDAO_H */
