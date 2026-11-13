#include <CGAL/SMDS_3/io_signature.h>

#include <cassert>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

int main()
{
  assert(CGAL::Get_io_signature<long>()() == "l");
  assert(CGAL::Get_io_signature<unsigned long>()() == "ul");
  assert(CGAL::Get_io_signature<long long>()() == "ll");
  assert(CGAL::Get_io_signature<unsigned long long>()() == "ull");

  std::string size_signature;
  if(std::is_same_v<std::size_t, unsigned int>)
    size_signature = "ui";
  else if(std::is_same_v<std::size_t, unsigned long>)
    size_signature = "ul";
  else
    size_signature = "ull";
  assert(CGAL::Get_io_signature<std::size_t>()() == size_signature);

  using Patch_index = std::pair<long, unsigned long>;
  using Index = std::variant<std::size_t, Patch_index, long long, unsigned long long>;
  assert(CGAL::Get_io_signature<Patch_index>()() == "std::pair<l,ul>");
  assert(CGAL::Get_io_signature<Index>()() ==
         "std::variant<" + size_signature + ",std::pair<l,ul>,ll,ull>");

  assert(CGAL::Get_io_signature<int>()() == "i");
  assert(CGAL::Get_io_signature<unsigned int>()() == "ui");
  assert(CGAL::Get_io_signature<short>()() == "s");
  assert(CGAL::Get_io_signature<unsigned short>()() == "us");
}
