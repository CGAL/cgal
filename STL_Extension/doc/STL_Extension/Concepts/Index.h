/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

The concept `Index` is a refinement of `Descriptor` which must be convertible from and to `std::size_t`, and must be incrementable and decrementable.

\cgalRefines `Descriptor`

\cgalHasModel int
\cgalHasModel size_t

\cgalHeading{Notation}

<dl>
<dt>`I`</dt>           <dd>Object of type Index.</dd>
<dt>`n`</dt>           <dd>Object of type `std::size_t`.</dd>
</dl>

\cgalHeading{Valid Expressions}

Expression                              | Returns                                                                      | Description
--------------------------------------- | ---------------------------------------------------------------------------- | -----------
`Index I(n)`                            | -                                                                            | Creates an index `I` from a `std::size_t`.
`(std::size_t)I`                        | `std::size_t`                                                                | Converts `I` into a `std::size_t`.
`++I`                                   | `Index&`                                                                     | Pre-increments `I`.
`I++`                                   | `Index`                                                                      | Post-increments `I`.
`--I`                                   | `Index&`                                                                     | Pre-decrements `I`.
`I--`                                   | `Index`                                                                      | Post-decrements `I`.
*/

class Index {};
