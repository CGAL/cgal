
#define LONG_NAME \
Wwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeee

template < class A >
struct LONG_NAME
{
  LONG_NAME (int i) : a(i) {}
  A a;
};

int main ()
{
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME<
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< int > > > > > > > > > > a (1);
  return 0;
}
