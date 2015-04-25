
template <typename FT, typename FT2>
void CHECK_EQUAL(const FT& a, const FT2& b)
{
  assert(a == b);
}

template <typename FT, typename FT2>
void CHECK_CLOSE(const FT& a, const FT2& b, const FT& eps)
{
  assert(CGAL::abs(a-b) < eps);
}
