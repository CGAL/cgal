
#ifndef TOOLS_H
#define TOOLS_H

#include <string>

#include <QVector2D>
#include <QVector3D>


std::string read_file(const std::string& file_name);

std::ostream& operator << (std::ostream& os, const QVector2D& v);
std::ostream& operator << (std::ostream& os, const QVector3D& v);
std::ostream& operator << (std::ostream& os, const QVector4D& v);


#endif
