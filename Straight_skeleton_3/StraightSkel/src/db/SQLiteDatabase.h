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
 * @file   db/SQLiteDatabase.h
 * @author Gernot Walzl
 * @date   2012-01-25
 */

#ifndef DB_SQLITEDATABASE_H
#define DB_SQLITEDATABASE_H

#include "db/ptrs.h"
#include <sqlite3.h>
#include <string>

namespace db {

class SQLiteDatabase {
public:
    SQLiteDatabase();
    virtual ~SQLiteDatabase();

    bool open(const std::string& filename);
    bool close();
    bool isOpened();

    SQLiteStmtSPtr prepare(const std::string& sql_query);

    bool beginTransaction();
    bool endTransaction();

protected:
    void printError();

    sqlite3* db_;
    bool transaction_;
};

}

#endif /* DB_SQLITEDATABASE_H */
