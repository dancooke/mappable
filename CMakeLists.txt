cmake_minimum_required(VERSION 3.5)
project(mappable)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MAPPABLE_SOURCES
    mappable/comparable.hpp
    mappable/contig_region.hpp
    mappable/genomic_region.hpp
    mappable/type_tricks.hpp
    mappable/mappable.hpp
    mappable/mappable_reference_wrapper.hpp
    mappable/mappable_range.hpp
    mappable/mappable_algorithms.hpp
    mappable/mappable_flat_multi_set.hpp
    mappable/mappable_flat_set.hpp
    mappable/mappable_map.hpp
    mappable/mappable_fwd.hpp
)

add_executable(example example.cpp ${MAPPABLE_SOURCES})
find_package (Boost 1.58)
if (Boost_FOUND)
    target_include_directories (example PUBLIC ${Boost_INCLUDE_DIR})
    target_link_libraries (example ${Boost_LIBRARIES})
endif (Boost_FOUND)
