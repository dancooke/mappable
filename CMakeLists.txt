cmake_minimum_required(VERSION 3.5)

project(mappable)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(MAPPABLE_SOURCES
    ${mappable_SOURCE_DIR}/mappable/comparable.hpp
    ${mappable_SOURCE_DIR}/mappable/contig_region.hpp
    ${mappable_SOURCE_DIR}/mappable/genomic_region.hpp
    ${mappable_SOURCE_DIR}/mappable/type_tricks.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_reference_wrapper.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_range.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_algorithms.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_flat_multi_set.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_flat_set.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_map.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_fwd.hpp
    ${mappable_SOURCE_DIR}/mappable/mappable_range_io.hpp
)

if (BUILD_TESTING)
    add_library(Mappable INTERFACE)
    target_sources(Mappable INTERFACE ${MAPPABLE_SOURCES})
    target_include_directories(Mappable INTERFACE ${mappable_SOURCE_DIR}/mappable)
    find_package (Boost 1.58)
    if (Boost_FOUND)
        target_include_directories (Mappable INTERFACE ${Boost_INCLUDE_DIR})
        target_link_libraries (Mappable ${Boost_LIBRARIES})
    endif (Boost_FOUND)
    add_subdirectory(test)
else()
    add_executable(example example.cpp ${MAPPABLE_SOURCES})
    find_package (Boost 1.58)
    if (Boost_FOUND)
        target_include_directories (example PUBLIC ${Boost_INCLUDE_DIR})
        target_link_libraries (example ${Boost_LIBRARIES})
    endif (Boost_FOUND)
    install(TARGETS example DESTINATION ${mappable_SOURCE_DIR})
endif(BUILD_TESTING)
