enable_testing()

# Generates a single test runner usable by CTest from a single test source file
function(add_boost_test SOURCE_FILE_NAME DEPENDENCY_LIB)
    # get the relative path in the source tree
    get_filename_component(test_path ${SOURCE_FILE_NAME} PATH)
    
    # get the source name without extension
    get_filename_component(test_name ${SOURCE_FILE_NAME} NAME_WE)
    
    # concatenate the relative path and name in an . separated identifier
    string(REPLACE "/" "." test_concat "${test_path}/${test_name}")
    
    # strip the trailing "_tests" part from the test ID
    string(REGEX REPLACE "_tests" "" test_id ${test_concat})
    
    # add a main for the tests
    add_definitions(-DBOOST_TEST_MAIN)
    
    add_executable(${test_id} ${SOURCE_FILE_NAME})

    target_link_libraries(${test_id} ${DEPENDENCY_LIB} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    
    # collect all the tests in the source file
    file(READ "${SOURCE_FILE_NAME}" SOURCE_FILE_CONTENTS)
    string(REGEX MATCHALL "BOOST_AUTO_TEST_CASE\\( *([A-Za-z_0-9]+) *\\)"
           FOUND_TESTS ${SOURCE_FILE_CONTENTS})
    
    foreach(HIT ${FOUND_TESTS})
        string(REGEX REPLACE ".*\\( *([A-Za-z_0-9]+) *\\).*" "\\1" unit_test_name ${HIT})
        
        add_test(NAME "${test_id}.${unit_test_name}"
                 COMMAND ${test_id}
                 --run_test=/*/*/${unit_test_name} --catch_system_error=yes)
    endforeach()
endfunction()

set(MAPPABLE_TEST_SOURCES
    unit_test_main.cpp
    comparable_tests.cpp
    contig_region_tests.cpp
    genomic_region_tests.cpp
    mappable_algorithm_tests.cpp
    mappable_flat_set_tests.cpp
    mappable_range_tests.cpp
    mappable_tests.cpp
)

add_definitions(-DBOOST_TEST_DYN_LINK)
find_package(Boost 1.58 REQUIRED COMPONENTS unit_test_framework REQUIRED)

include_directories(${Boost_INCLUDE_DIRS} ${mappable_SOURCE_DIR}/mappable ${mappable_SOURCE_DIR}/test)

set(TEST_DEPENDENCY_LIBS Mappable)

# Add each test
foreach(SRC ${MAPPABLE_TEST_SOURCES})
    add_boost_test(${SRC} "${TEST_DEPENDENCY_LIBS}")
endforeach()
