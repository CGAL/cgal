namespace CGAL {

/*!
\ingroup PkgMatrixSearch

\brief returns the element `x` in one of the sorted matrices from the 
range `[f, l)`, for which `t.is_feasible(x)` 
is true and `t.compare(x, y)` is true for all other 
`y` values from any matrix for which `t.is_feasible(y)` is true. 


The function `sorted_matrix_search()` selects the smallest entry 
in a set of sorted matrices that fulfills a certain feasibility 
criterion. 

More exactly, a matrix \f$ M = (m_{i j}) \in S^{r \times l}\f$ 
(over a totally ordered set \f$ S\f$) is sorted, iff 
\f{eqnarray*}{
\forall \, 1 \le i \le r,\; 1 \le j < l\; :\; m_{i j} \le m_{i (j+1)} 
\;\; {\it and}\\ 
\forall \, 1 \le i < r,\; 1 \le j \le l\; :\; m_{i j} \le m_{(i+1) j} 
\;\;. 
\f}

Now let \f$ \mathcal{M}\f$ be a set of \f$ n\f$ sorted matrices over \f$ S\f$ 
and \f$ f\f$ be a monotone predicate on \f$ S\f$, i.e.\ 
\f[ 
f\: :\: S \longrightarrow\, \textit{bool} \quad{\rm with}\quad f(r) 
\;\Longrightarrow\; \forall\, t \in S\,,\: t > r \; :\; f(t)\;. 
\f] 

If we assume there is any feasible element in one of the matrices 
in \f$ \mathcal{M}\f$, there certainly is a smallest such element. This is the one 
we are searching for. 

The feasibility test as well as some other parameters can (and 
have to) be customized through a traits class. 


\pre <OL> 
<LI>All matrices in \f$ \left[f,\, l\right)\f$ are sorted according 
to `Traits::compare_non_strictly`. 
<LI>There is at least one entry \f$ x\f$ in a matrix \f$ M \in 
\left[f,\, l\right)\f$ for which `Traits::is_feasible(x)` is 
true. 
</OL> 

\tparam Traits is a model for `SortedMatrixSearchTraits`.
\tparam RandomAccessIterator has `Traits::Matrix` as value type.

\cgalHeading{Implementation}

The implementation uses an algorithm by 
Frederickson and Johnson\cgalCite{fj-fkppc-83}, \cgalCite{fj-gsrsm-84} and runs in 
\f$ \mathcal{O}(n \cdot k + f \cdot \log (n \cdot k))\f$, where \f$ n\f$ is 
the number of input matrices, \f$ k\f$ denotes the maximal dimension of 
any input matrix and \f$ f\f$ the time needed for one feasibility test. 


\sa `SortedMatrixSearchTraits` 
*/
template < class RandomAccessIterator, class
Traits > Traits::Value sorted_matrix_search(
RandomAccessIterator f, RandomAccessIterator l, const Traits&
t);

} /* namespace CGAL */

