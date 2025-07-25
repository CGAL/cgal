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
 * @file   db/SQLiteDatabase.cpp
 * @author Gernot Walzl
 * @date   2012-01-25
 */

#include "db/SQLiteDatabase.h"

#include "debug.h"
#include "db/SQLiteStmt.h"
#include "util/StringFactory.h"
#include <iostream>

namespace db {

SQLiteDatabase::SQLiteDatabase() {
    this->db_ = nullptr;
    this->transaction_ = false;
}

SQLiteDatabase::~SQLiteDatabase() {
    this->close();
}

bool SQLiteDatabase::open(const std::string& filename) {
    CGAL_SS3_IO_TRACE(filename);
    bool result = false;
    if (SQLITE_OK == sqlite3_open(filename.c_str(), &db_)) {
        result = true;
    } else {
        this->printError();
        this->db_ = nullptr;
    }
    return result;
}

bool SQLiteDatabase::close() {
    bool result = false;
    if (db_) {
        if (SQLITE_OK == sqlite3_close(db_)) {
            db_ = nullptr;
            result = true;
        }
    }
    return result;
}

bool SQLiteDatabase::isOpened() {
    return (db_ != nullptr);
}

void SQLiteDatabase::printError() {
    if (db_) {
        std::cout << sqlite3_errmsg(db_) << std::endl;
    }
}

SQLiteStmtSPtr SQLiteDatabase::prepare(const std::string& sql_query) {
    SQLiteStmtSPtr result = SQLiteStmtSPtr();
    if (db_) {
        sqlite3_stmt* stmt;
        CGAL_SS3_IO_TRACE(sql_query);
        if (SQLITE_OK == sqlite3_prepare(db_, sql_query.c_str(), -1, &stmt, NULL)) {
            result = SQLiteStmtSPtr(new SQLiteStmt(db_, stmt));
        } else {
            this->printError();
        }
    }
    return result;
}

bool SQLiteDatabase::beginTransaction() {
    bool result = false;
    if (!transaction_) {
        CGAL_SS3_IO_TRACE("BEGIN TRANSACTION;");
        if (SQLITE_OK == sqlite3_exec(db_, "BEGIN TRANSACTION;",
                NULL, NULL, NULL)) {
            result = true;
            transaction_ = true;
        } else {
            this->printError();
        }
    }
    return result;
}

bool SQLiteDatabase::endTransaction() {
    bool result = false;
    if (transaction_) {
        CGAL_SS3_IO_TRACE("END TRANSACTION;");
        if (SQLITE_OK == sqlite3_exec(db_, "END TRANSACTION;",
                NULL, NULL, NULL)) {
            result = true;
            transaction_ = false;
        } else {
            this->printError();
        }
    }
    return result;
}

}
