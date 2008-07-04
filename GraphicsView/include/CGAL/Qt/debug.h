#ifndef CGAL_QT_DEBUG_H
#define CGAL_QT_DEBUG_H

#include <QString>

namespace CGAL {
namespace Qt {

/**
 *  Must be used like that:
 *     CGAL::Qt:traverse_resources(":/cgal"); // view CGAL resources
 *  or
 *     CGAL::Qt:traverse_resources(":"); // view all resources
 *  and displays the resources tree on std::cerr.
 */
void traverse_resources(const QString& name,
                        const QString& dirname = QString(),
                        int indent = 0);

} // namespace Qt
} // namespace CGAL


#endif // CGAL_QT_DEBUG_H
