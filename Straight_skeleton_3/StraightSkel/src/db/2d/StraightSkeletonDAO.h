/**
 * @file   db/2d/StraightSkeletonDAO.h
 * @author Gernot Walzl
 * @date   2013-05-23
 */

#ifndef DB_2D_STRAIGHTSKELETONDAO_H
#define DB_2D_STRAIGHTSKELETONDAO_H

#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/StraightSkeleton.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::skel::StraightSkeleton;
using data::_2d::skel::StraightSkeletonSPtr;

class StraightSkeletonDAO {
friend class DAOFactory;
public:
    virtual ~StraightSkeletonDAO();
    std::string getTableSchema() const;
    int createSkelID(StraightSkeletonSPtr skel);
    int insert(StraightSkeletonSPtr skel);
    bool del(StraightSkeletonSPtr skel);
    StraightSkeletonSPtr find(int skelid);
    int findPolyID(int skelid);
    bool update(StraightSkeletonSPtr skel);
private:
    StraightSkeletonDAO();
    int nextSkelID();
};

} }

#endif /* DB_2D_STRAIGHTSKELETONDAO_H */
