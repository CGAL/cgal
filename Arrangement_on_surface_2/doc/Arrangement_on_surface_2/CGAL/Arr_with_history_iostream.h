namespace CGAL {

/*!
\ingroup PkgArrangement2Read

Reads the arrangement-with-history object `arr` from the given 
input stream `is` using a specific input format defined by
\"formatter\".

*/
template<typename Traits, typename Dcel,
         typename WithHistoryFormatter>
std::istream& read (Arrangement_with_history_2<Traits,Dcel>& arr,
                    std::istream& is,
                    WithHistoryFormatter& formatter);

/*!
\ingroup PkgArrangement2Write
Writes the arrangement-with-history object `arr` into the given
output stream `os` using a specific output format defined by
`formatter`.
*/
template<typename Traits, typename Dcel,
         typename WithHistoryFormatter>
std::ostream& write (const Arrangement_with_history_2<Traits,Dcel>& arr,
                     std::ostream& os,
                     WithHistoryFormatter& formatter);
 
/*!
\ingroup PkgArrangement2op_left_shift
Inserts the arrangement-with-history object `arr` into the output
stream `os` using the output format defined by the
`Arr_with_history_text_formatter` class. Only the basic geometric
and topological features of the arrangement are inserted. Auxiliary
data that may be attached to the <span class="textsc">Dcel</span> features is ignored.
*/
template<typename Traits, typename Dcel>
std::ostream& operator<< (std::ostream& os,
                          const Arrangement_with_history_2<Traits,Dcel>& arr);

/*!
\ingroup PkgArrangement2op_right_shift
Extracts an arrangement-with-history from a given input stream using
the default input format.
*/
template<class Traits, class Dcel>
std::istream& operator>>(std::istream& is, Arrangement_with_history_2<Traits,Dcel>& arr);
}
