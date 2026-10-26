/*
Copyright (C) 2026  TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <vector>

namespace CGAL::internal {

template <class QT>
concept QueueTraits = requires(QT traits, typename QT::Element_handle elt, int i) {
	typename QT::Element_handle;

	{traits.setIndex(elt, i)};

	{traits.getIndex(elt)}->std::same_as<int>;

	{traits.compare(elt, elt)}
	    ->std::same_as<int>; // negative if elt < elt2, positive if elt > elt2, zero if elt = elt2. Smallest value == highest priority (top of queue)
};

template <QueueTraits QT> class IndexedPriorityQueue {
  public:
	using Element_handle = QT::Element_handle;

  private:
	std::vector<Element_handle> queue;
	QT traits;

	void siftUp(int k, Element_handle elt);
	void siftDown(int k, Element_handle elt);

  public:
	IndexedPriorityQueue() requires std::is_default_constructible<QT>::value {}
	IndexedPriorityQueue(QT traits) : traits(std::move(traits)) {}

	bool empty() const;

	void push(Element_handle elt);
	Element_handle pop();
	Element_handle peek() const;

	bool remove(Element_handle elt);
	bool contains(Element_handle elt) const;
	void update(Element_handle elt);

	void clear();

	const std::vector<Element_handle>& content() const;
};

} // namespace CGAL::internal

#include "indexed_priority_queue.hpp"