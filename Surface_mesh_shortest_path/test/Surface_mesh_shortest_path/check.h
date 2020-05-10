
template <typename FT, typename FT2>
void CHECK_EQUAL(const FT& a, const FT2& b)
{
  if (a != b)
    std::cerr << "ERROR: a (" << a << ") is not equal to b (" << b << ").\n";

  assert(a == b);
}

template <typename FT, typename FT2>
void CHECK_CLOSE(const FT& a, const FT2& b, const FT& eps)
{
  if ( !(CGAL::abs(a-b) < eps) )
    std::cerr << "ERROR: difference (" << CGAL::abs(a-b) << ") is larger than eps (" << eps << ").\n";
  assert(CGAL::abs(a-b) < eps);
}
