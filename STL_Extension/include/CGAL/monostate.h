#ifndef CGAL_MONOSTATE_H
#define CGAL_MONOSTATE_H

namespace CGAL {

struct monostate {
  void* ptr = nullptr;
};

constexpr bool operator==(monostate, monostate) noexcept
{
  return true;
}

constexpr bool operator!=(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator<(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator>(monostate, monostate) noexcept
{
  return false;
}

constexpr bool operator<=(monostate, monostate) noexcept
{
  return true;
}

constexpr bool operator>=(monostate, monostate) noexcept
{
  return true;
}

} // namespace CGAL

#endif