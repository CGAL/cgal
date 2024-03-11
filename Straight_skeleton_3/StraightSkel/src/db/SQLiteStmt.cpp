/**
 * @file   db/SQLiteStmt.cpp
 * @author Gernot Walzl
 * @date   2012-01-31
 */

#include "db/SQLiteStmt.h"

#include <iostream>

namespace db {

SQLiteStmt::SQLiteStmt(sqlite3* db, sqlite3_stmt* stmt) {
    this->db_ = db;
    this->stmt_ = stmt;
    this->result_set_ = false;
    this->result_set_row_ = 0;
}

SQLiteStmt::~SQLiteStmt() {
    this->close();
}

bool SQLiteStmt::close() {
    bool result = false;
    if (stmt_) {
        if (SQLITE_OK == sqlite3_finalize(stmt_)) {
            result = true;
            stmt_ = nullptr;
            result_set_ = false;
            result_set_row_ = 0;
        }
    }
    return result;
}

void SQLiteStmt::printError() {
    if (db_) {
        std::cout << sqlite3_errmsg(db_) << std::endl;
    }
}

bool SQLiteStmt::bindInteger(int col, int value) {
    bool result = false;
    if (db_ && stmt_) {
        if (SQLITE_OK == sqlite3_bind_int(stmt_, col, value)) {
            result = true;
        } else {
            this->printError();
        }
    }
    return result;
}

bool SQLiteStmt::bindDouble(int col, double value) {
    bool result = false;
    if (db_ && stmt_) {
        if (SQLITE_OK == sqlite3_bind_double(stmt_, col, value)) {
            result = true;
        } else {
            this->printError();
        }
    }
    return result;
}

bool SQLiteStmt::bindString(int col, std::string value) {
    bool result = false;
    if (db_ && stmt_) {
        if (SQLITE_OK == sqlite3_bind_text(stmt_, col, value.c_str(), -1, SQLITE_STATIC)) {
            result = true;
        } else {
            this->printError();
        }
    }
    return result;
}

int SQLiteStmt::execute() {
    int result = 0;
    result_set_ = false;
    if (db_ && stmt_) {
        int returned = sqlite3_step(stmt_);
        if (returned == SQLITE_ERROR) {
            this->printError();
        } else {
            result = 1;
            if (returned == SQLITE_ROW) {
                result_set_ = true;
            } else {
                result = sqlite3_changes(db_);
            }
        }
    }
    return result;
}

bool SQLiteStmt::fetchRow() {
    bool result = false;
    if (!result_set_) {
        return false;
    }
    if (this->result_set_row_ == 0) {
        this->result_set_row_++;
        result = true;
    } else if (this->result_set_row_ > 0) {
        if (SQLITE_ROW == sqlite3_step(stmt_)) {
            this->result_set_row_++;
            result = true;
        }
    }
    return result;
}

int SQLiteStmt::getInteger(int col) {
    return sqlite3_column_int(stmt_, col);
}

double SQLiteStmt::getDouble(int col) {
    return sqlite3_column_double(stmt_, col);
}

std::string SQLiteStmt::getString(int col) {
    return std::string((const char*)sqlite3_column_text(stmt_, col));
}

bool SQLiteStmt::reset() {
    bool result = false;
    result = (sqlite3_reset(stmt_) == SQLITE_OK);
    return result;
}

}
