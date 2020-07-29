#ifndef MYPOINTC2_IOSTREAM_H
#define MYPOINTC2_IOSTREAM_H

std::ostream &
operator<<(std::ostream &os, const MyPointC2 &p)
{
    switch(CGAL::get_mode(os)) {
    case CGAL::IO::ASCII :
        return os << p.x() << ' ' << p.y() << ' ' << p.color();
    case CGAL::IO::BINARY :
        CGAL::write(os, p.x());
        CGAL::write(os, p.y());
        CGAL::write(os, p.color());
        return os;
    default:
        return os << "MyPointC2(" << p.x() << ", " << p.y() << ", " << p.color() << ')';
    }
}



std::istream &
operator>>(std::istream &is, MyPointC2 &p)
{
    double x, y;
    int c;
    switch(CGAL::get_mode(is)) {
    case CGAL::IO::ASCII :
      is >> x >> y >> c;
        break;
    case CGAL::IO::BINARY :
        CGAL::read(is, x);
        CGAL::read(is, y);
        CGAL::read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is) {
      p = MyPointC2(x, y, c);
    }
    return is;
}
#endif //MYPOINTC2_IOSTREAM_H
