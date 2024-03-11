/**
 * @file   db/SQLiteStmt.h
 * @author Gernot Walzl
 * @date   2012-01-31
 */

#ifndef DB_SQLITESTMT_H
#define DB_SQLITESTMT_H

#include <sqlite3.h>
#include <string>

namespace db {

class SQLiteStmt {
public:
    SQLiteStmt(sqlite3* db, sqlite3_stmt* stmt);
    virtual ~SQLiteStmt();

    bool close();

    /**
     * The leftmost SQL parameter has an index of 1.
     */
    bool bindInteger(int col, int value);
    bool bindDouble(int col, double value);
    bool bindString(int col, std::string value);
    int execute();

    bool fetchRow();
    /**
     * The leftmost column of the result set has the index 0.
     */
    int getInteger(int col);
    double getDouble(int col);
    std::string getString(int col);

    bool reset();

protected:
    void printError();
    sqlite3* db_;
    sqlite3_stmt* stmt_;
    bool result_set_;
    unsigned int result_set_row_;
};

}

#endif /* DB_SQLITESTMT_H */
