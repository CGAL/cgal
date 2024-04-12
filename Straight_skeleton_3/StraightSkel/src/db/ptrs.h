// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

/**
 * @file   db/ptrs.h
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#ifndef DB_PTRS_H
#define DB_PTRS_H

#include "smarter_ptr.h"

namespace db {

class SQLiteDatabase;
class SQLiteStmt;

typedef SHARED_PTR<SQLiteDatabase> SQLiteDatabaseSPtr;
typedef WEAK_PTR<SQLiteDatabase> SQLiteDatabaseWPtr;
typedef SHARED_PTR<SQLiteStmt> SQLiteStmtSPtr;
typedef WEAK_PTR<SQLiteStmt> SQLiteStmtWPtr;

}

#endif /* DB_PTRS_H */
