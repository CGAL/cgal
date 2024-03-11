/**
 * @file   sqlite_sample.c
 * @author Gernot Walzl
 * @date   2012-01-25
 */

/* compiling and linking
 gcc -lsqlite3 -ldl sqlite_sample.c
 */

/* creating data
sqlite3 test.db3 "CREATE TABLE mytable (firstcolumn TEXT);"
sqlite3 test.db3 "INSERT INTO mytable VALUES ('This is a simple test.');"
sqlite3 test.db3 "SELECT * FROM mytable;"
 */

#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>

int main(int argc, char** argv) {
    sqlite3* db;
    if (SQLITE_OK != sqlite3_open("test.db3", &db)) {
        printf("Error on open.\n");
        return (EXIT_FAILURE);
    }
    sqlite3_stmt* stmt;
    if (SQLITE_OK != sqlite3_prepare(db, "SELECT * FROM mytable", -1, &stmt, NULL)) {
        printf("Error on prepare.\n");
        return (EXIT_FAILURE);
    }
    while (SQLITE_ROW == sqlite3_step(stmt)) {
        const unsigned char * value = sqlite3_column_text(stmt, 0);
        printf("%s\n", value);
    }
    if (SQLITE_OK == sqlite3_finalize(stmt)) {
        sqlite3_close(db);
        printf("Exit gracefully.\n");
    }
    return (EXIT_SUCCESS);
}
