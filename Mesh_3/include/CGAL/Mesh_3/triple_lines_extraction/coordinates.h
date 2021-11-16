#ifndef COORDINATES_H
#define COORDINATES_H

#include <array>

typedef std::array<int, 3> Coordinates;
constexpr Coordinates coordinates[8] = { 0, 0, 0,
                                         1, 0, 0,
                                         0, 1, 0,
                                         1, 1, 0,
                                         0, 0, 1,
                                         1, 0, 1,
                                         0, 1, 1,
                                         1, 1, 1 };

Coordinates operator-(Coordinates b, Coordinates a) {
  return { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
}

Coordinates cross(Coordinates a, Coordinates b) {
  return { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
}

Coordinates square(Coordinates c) {
  return { c[0]*c[0], c[1]*c[1], c[2]*c[2] };
}

int dist(Coordinates a, Coordinates b) {
  auto s = square(b - a);
  return s[0] + s[1] + s[2];
}

#endif // COORDINATES_H
