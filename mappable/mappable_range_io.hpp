// Copyright (c) 2017-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mappable_range_io_hpp
#define mappable_range_io_hpp

#include <string>
#include <iterator>
#include <algorithm>
#include <iostream>

#include "mappable_range.hpp"

namespace mappable {

namespace detail {

template <typename T, typename Range>
std::ostream& pretty_print_range(std::ostream& os, const Range& range, const char* delim)
{
    if (!empty(range)) {
        std::copy(std::cbegin(range), std::prev(std::cend(range)), std::ostream_iterator<T> {os, delim});
        os << range.back();
    }
    return os;
}

} // namespace detail

static std::string range_io_delim = " ";

template <typename It>
std::ostream& operator<<(std::ostream& os, const OverlapRange<It>& range)
{
    using T = typename std::iterator_traits<It>::value_type;
    return detail::pretty_print_range<T>(os, range, range_io_delim.c_str());
}

template <typename It>
std::ostream& operator<<(std::ostream& os, const ContainedRange<It>& range)
{
    using T = typename std::iterator_traits<It>::value_type;
    return detail::pretty_print_range<T>(os, range, range_io_delim.c_str());
}

} // namespace mappable

#endif
