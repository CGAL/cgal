
struct Class {
  struct Nested_class {
    Nested_class(){}
  };
};


template <class Cl>
void
function(const Cl& cl,
	 typename Cl::Nested_class start =  typename Cl::Nested_class())
{}


int
main()
{
  Class cl;
  function(cl);
  return 0;
}
