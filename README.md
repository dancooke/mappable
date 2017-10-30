# Mappable

Mappable is a powerful C++ template library for genomic region storage, manipulation, and querying. Mappable is independent of the data structure used to store the underlying objects, instead relying on C++'s [curiously recurring template pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) to define region concepts and ordering. This allows Mappable to build on standard C++ containers and algorithms to give an efficient, flexible, and expressive set of region based containers and algorithms which act directly on `Mappable` objects. Mappable has been heavily tested and benchmarked, and can perform orders of magnitude faster than other approaches (see `benchmark.cpp`).

## Requirements

* A C++14 compiler and standard library implementation
* Boost 1.58 or greater

To compile the example you will also need CMake 3.5 or greater:

```shell
$ cmake .
$ make
$ ./example
```

## Basic usage

The starting point for any `Mappable` type is to inherit from `Mappable` and implement the `mapped_region` member method:

    struct MyMappableType : public Mappable<MyMappableType>
    {
        ContigRegion region;
        MyDataType data;
        const ContigRegion& mapped_region();
        MyMappableType(ContigRegion region, MyDataType data);
    };

Now `MyMappableType` can be used with all of Mappable's algorithms:

    std::vector<MyMappableType> mappables {};
    mappables.emplace_back(ContigRegion {0, 1}, "A");
    mappables.emplace_back(ContigRegion {1, 3}, "B");
    mappables.emplace_back(ContigRegion {2, 4}, "C");
    std::cout << count_contained(mappables, ContigRegion {1, 3}) << std::endl; // 2
    std::cout << count_overlapped(mappables, ContigRegion {1, 3}) << std::endl; // 1

The library also includes some helpful containers that are heavily optimised for `Mappable` types:

    MyMappableType a {ContigRegion {0, 1}, "A"}, b {ContigRegion {1, 3}, "B"};
    MappableFlatSet<MyMappableType> mappables {a};
    mappables.insert(b);
    mappables.emplace(ContigRegion {1, 3}, "C");
    mappables.erase(b);
