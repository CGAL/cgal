namespace A {

class Cls {};

bool
foo(const Cls&)
{ return true; }

} // namespace A

bool
foo(int)
{ return false; }

int
main()
{
  A::Cls c;
  foo( c );   // Koenig lookup finds A::foo(const Cls&)
  return 0;
}
