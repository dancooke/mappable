# Mappable

Mappable is a powerful C++ template library for genomic region storage, manipulation, and querying. Mappable is independent of the data structure used to store the underlying objects, instead relying on C++'s [curiously recurring template pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to define region concepts and ordering. This allows Mappable to build on standard C++ containers and algorithms to give an efficient, flexible, and expressive set of region based containers and algorithms which act directly on `Mappable` objects. Mappable was originally developed as part of the variant caller [octopus](https://github.com/luntergroup/octopus), where it is used extensively. It has been heavily tested and benchmarked, and can perform orders of magnitude faster than other approaches (see `benchmark.cpp`).

## Requirements

* A C++14 compiler and standard library implementation
* Boost 1.58 or greater

To compile the example you will also need CMake 3.5 or greater:

```shell
$ mkdir build && cd build
$ cmake .. && make
$ ./example
```

## Basic usage

The starting point for any `Mappable` type is to inherit from `Mappable` and implement the `mapped_region` member method:

```cpp
struct MyMappableType : public Mappable<MyMappableType>
{
    ContigRegion region;
    MyDataType data;
    const auto& mapped_region() const noexcept { return region; }; // required by Mappable
    // Constructors and other member methods
};
```

Now `MyMappableType` can be used with all of Mappable's algorithms:

```cpp
std::vector<MyMappableType> mappables {};
mappables.emplace_back(ContigRegion {0, 1}, "A");
mappables.emplace_back(ContigRegion {1, 3}, "B");
mappables.emplace_back(ContigRegion {2, 4}, "C");
auto overlapped = overlap_range(mappables, ContigRegion {1, 3}); // iterator range
std::cout << overlapped << std::endl;
std::cout << count_contained(mappables, ContigRegion {1, 3}) << std::endl; // 1
```

The library also includes some helpful containers that are heavily optimised for `Mappable` types:

```cpp
MyMappableType a {ContigRegion {0, 1}, "A"}, b {ContigRegion {1, 3}, "B"};
MappableFlatSet<MyMappableType> mappables {a};
mappables.insert(b);
mappables.emplace(ContigRegion {1, 3}, "C");
mappables.erase(b);
```

## How it works

Fundamental to a `Mappable` object is the `mapped_region()` public method which enables compile time static polymorphism using the CRTP. As the base `Mappable` template class is empty, there is no memory overhead to deriving from `Mappable`. The `mapped_region()` method must return one of the two fundamental region types, `ContigRegion` or `GenomicRegion`. `ContigRegion` simply defines a range, whilst `GenomicRegion` defines a contig plus `ContigRegion`. Both `ContigRegion` and `GenomicRegion` have a comprehensive free function interface which can be invoked by any `Mappable` object via a callback to `mapped_region()`. For example, both  `ContigRegion` and `GenomicRegion` implement the `overlaps(const M&lhs, const M& rhs)` method which returns whether or not `lhs` and `rhs` are overlapping. In `Mappable`s interface, there is also an `overlaps` method that accepts any two `Mappable` types that forwards the call to the correct `overlaps` implementation. The first benefit of this approach is it allows users to call `overlaps` using `overlaps(lhs, rhs)` rather than `overlaps(lhs.mapped_region(), rhs.mapped_region())`. However, more importantly, it allows one to write general algorithms that are independent of the `Mappable` objects they operate on.

The Mappable algorithms library is designed in a similar way to the C++ STL; most algorithms operate on collections of `Mappable` objects via iterator ranges. Usually the range must be sorted with respect to `operator<` for the underlying `Mappable` region type (`ContigRegion` or `GenomicRegion`). How the `Mappable` objects are actually stored in memory is only important so far as it determines the [iterator type](http://en.cppreference.com/w/cpp/concept/Iterator) that each algorithm can operate on. This allows users to choose a data structure that has the best time and memory complexity requirements for their particular use case.

Although mappable algorithms can be used directly with standard conforming containers, the Mappable library includes some of it's own containers. The reason for this is that some mappable algorithms (predominantly `overlap_range`) have faster implementations if the range of `Mappable` objects satisfies a stronger ordering than that defined by `operator<`. In particular, `operator<`implies an ordering from 'left to right', but does not, and in the general case cannot, require a 'right to left' ordering. However, some collections will naturally satisfy the 'right to left' ordering in addition to the 'left to right' one required by `operator<`. If this is the case,`overlap_range` has worst-case complexity logarithmic (apposed to amortized worse case logarithmic). There is also a faster implementation of `overlap_range` if the maximum sized`Mappable` object is known. The provided containers simply keep track of these properties, and call the most appropriate mappable algorithm.
