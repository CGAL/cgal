/*!
\cgalConcept

A CombinationElement can be used as template parameter for the class
`Combination_enumerator<CombinationElement>`.

\cgalHasModelsBegin
\cgalHasModelsBare{Any integer type (`char`, `short`, `int`, `long`, etc.)}
\cgalHasModelsBare{Pointers}
\cgalHasModelsBare{Random access iterators}
\cgalHasModelsEnd

\sa `CGAL::Combination_enumerator<CombinationElement>`


*/

class CombinationElement {
public:

/// \name Creation
/// @{

/*!
Copy constructor
*/
CombinationElement(const CombinationElement & e2);

/// @}


/// \name Types
/// @{

/*!
the type of point being generated.
*/
typedef unspecified_type value_type;

/// @}

/// \name Operations
/// @{

/*!
Incrementation
*/
void operator++();

/*!
Decrementation
*/
void operator--();

/*!
Total order comparison
*/
bool operator<(const CombinationElement & e2);

/*!
Equality test
*/
bool operator==(const CombinationElement & e2);


/*!
Equivalent to calling `++(*this)` `i` times if i is positive.
Equivalent to calling `--(*this)` `-i` times if i is negative.
*/
CombinationElement operator+(int i);

/*!
Equivalent to calling ++(*this) `i` times if i is positive.
Equivalent to calling --(*this) `-i` times if i is negative.
*/
int operator-(const CombinationElement & e2);

/// @}

}; /* end CombinationElement */
