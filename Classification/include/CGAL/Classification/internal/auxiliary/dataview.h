// Copyright (c) 2014 Stefan Walk
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: MIT
//
// Author(s)     : Stefan Walk

// Modifications from original library:
//  * changed inclusion protection tag
//  * moved to namespace CGAL::internal::


#ifndef CGAL_INTERNAL_LIBLEARNING_DATAVIEW_H
#define CGAL_INTERNAL_LIBLEARNING_DATAVIEW_H
#include <vector>

namespace CGAL { namespace internal {
    
namespace liblearning {
//! \brief A view for one-dimensional data
//
// The data can be e.g. a row or column of a matrix or the elements of an
// array or a std::vector - anything with a fixed step size between pointers
template <typename ElementType>
struct DataView1D {
    //! \brief Element access
    //
    // \param idx the index of the element
    ElementType& operator()(size_t idx) const
    {
        return *(data + step * idx);
    }
    //! \brief Construct using pointer, size and (optional) step size
    DataView1D(ElementType* ptr, size_t size, ptrdiff_t step_size = 1) :
        data(ptr), step(step_size), num_elements(size)
    {
    }
    //! \brief Construct view from std::vector
    //
    // \param vec the vector, by-ref since we don't take ownership - no
    //            temporaries allowed. maybe needs to change
    DataView1D(std::vector<ElementType>& vec) :
        data(&vec[0]), step(1), num_elements(vec.size())
    {
    }
    
    //! \brief Construct empty view
    DataView1D() : data(0), step(1), num_elements(0)
    {
    }
    ElementType* data; //!< Pointer to first element
    ptrdiff_t step;    //!< Step size between elements - for a vector this is 1, no need to multiply by sizeof(ElementType)
    size_t num_elements; //!< Number of elements in the view
};

template <typename ElementType>
struct DataView2D {
    //! \brief Construct empty view
    DataView2D() : data(0), row_step(1), col_step(1), rows(0), cols(0)
    {}
    //! \brief Construct view from memory using given step sizes
    DataView2D(ElementType* ptr, size_t rows_, size_t cols_, ptrdiff_t row_step_, ptrdiff_t col_step_) :
        data(ptr), row_step(row_step_), col_step(col_step_), rows(rows_), cols(cols_)
    {
    }
    //! \brief Construct view from a continuous block of memory in row-major order
    DataView2D(ElementType* ptr, size_t rows_, size_t cols_) :
        data(ptr), row_step(cols_), col_step(1), rows(rows_), cols(cols_)
    {
    }
    //! \brief Element access
    ElementType& operator()(size_t row_idx, size_t col_idx) const
    {
        return *(data + row_step * row_idx + col_step * col_idx);
    }
    //! \brief Return a 1D view of a row
    DataView1D<ElementType> row(size_t row_idx)
    {
        return DataView1D<ElementType>(data + row_step * row_idx, cols, col_step);
    }
    //! \brief Return a 1D view of a column
    DataView1D<ElementType> col(size_t col_idx)
    {
        return DataView1D<ElementType>(data + col_step * col_idx, rows, row_step);
    }
    //! \brief Return a new view, using a subset of rows
    DataView2D row_range(size_t start_row, size_t end_row) const
    {
        DataView2D ret(*this);
        ret.data = data + row_step * start_row;
        ret.rows = end_row - start_row;
        return ret;
    }
    //! \brief Return a new view, using a subset of columns
    DataView2D col_range(size_t start_col, size_t end_col) const
    {
        DataView2D ret(*this);
        ret.data = data + col_step * start_col;
        ret.cols = end_col - start_col;
        return ret;
    }
    //! \brief Transpose the matrix (actually done by swapping the steps, no
    //                               copying)
    DataView2D transpose() const
    {
        DataView2D ret(*this);
        std::swap(ret.row_step, ret.col_step);
        std::swap(ret.rows,     ret.cols);
        return ret;
    }
    //! \brief Create a 2D view from a 1D one, represent as a column
    static DataView2D column(DataView1D<ElementType> vec)
    {
        return DataView2D(vec.data, vec.num_elements, 1);
    }
    //! \brief Create a 2D view from a 1D one, represent as a row
    static DataView2D row(DataView1D<ElementType> vec)
    {
        return DataView2D(vec.data, 1, vec.num_elements);
    }
    //! \brief Return the number of elements in this view
    size_t num_elements() { return rows * cols; }
    //! \brief Return true if the view is empty (no elements)
    bool   empty() { return num_elements() == 0; }
    //! \brief Return true if rows in this view are continuous in memory
    bool row_continuous() {
        return col_step == 1;
    }
    //! \brief Return true if columns in this view are continuous in memory
    bool col_continuous() {
        return row_step == 1;
    }
    //! \brief Get pointer to row (only valid if row-continuous)
    ElementType* row_pointer(size_t row_idx) {
        return data + row_idx * row_step;
    }
    //! \brief Get pointer to column (only valid if column-continuous)
    ElementType* col_pointer(size_t col_idx) {
        return data + col_idx * col_step;
    }

    ElementType* data; //!< Pointer to first element
    ptrdiff_t row_step; //!< Pointer difference between an element and its right neighbor - 1 if row-continuous
    ptrdiff_t col_step; //!< Pointer difference between an element and its bottom neighbor - 1 if column-continuous
    size_t rows; //!< Number of rows in the view
    size_t cols; //!< Number of columns in the view
};

//! \brief Determine if two 2D views have the same dimensions
template <typename A, typename B>
bool equal_dims(DataView2D<A> view_a, DataView2D<B> view_b)
{
    return view_a.rows == view_b.rows && view_a.cols == view_b.cols;
}
}

}} // namespace CGAL::internal::

#endif
