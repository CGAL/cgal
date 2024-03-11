/**
 * @file   db/2d/ArcDAO.h
 * @author Gernot Walzl
 * @date   2013-05-23
 */

#ifndef DB_2D_ARCDAO_H
#define DB_2D_ARCDAO_H

#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/Arc.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::skel::Arc;
using data::_2d::skel::ArcSPtr;

class ArcDAO {
friend class DAOFactory;
public:
    virtual ~ArcDAO();
    std::string getTableSchema() const;
    int insert(ArcSPtr arc);
    bool del(ArcSPtr arc);
    ArcSPtr find(int skelid, int aid);
    bool update(ArcSPtr arc);
protected:
    ArcDAO();
    int nextAID(int skelid);
};

} }

#endif /* DB_2D_ARCDAO_H */

