#ifndef FLAVORED_OBJECT_H
#define FLAVORED_OBJECT_H

#include <CEP/Vanilla/flavor.h>

template <class Object>
class Flavored_object : public Object
{
public:
   Flavored_object()
   { }

   Flavored_object(Flavor f) : _flavor(f)
   { }

   Flavored_object(Object o) : Object(o), _flavor(VANILLA)
   { }

   Flavored_object(Object o, Flavor f): Object(o), _flavor(f)
   { }

   void set_flavor(Flavor f)
   {
      _flavor = f;
   }

   void enhance_flavor()
   {
      _flavor = flavor_enhance(_flavor);
   }

   bool is_valid() const
   {
      return valid_flavor(_flavor);
   }

   Flavor flavor() const
   {
      return _flavor;
   }

private:
   Flavor _flavor;
};

#endif // FLAVORED_OBJECT_H
