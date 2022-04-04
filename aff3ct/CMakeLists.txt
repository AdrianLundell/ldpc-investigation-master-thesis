cmake_minimum_required(VERSION 3.13) # Need FetchContent
cmake_policy(SET CMP0054 NEW)

# Compilation options
option(AFF3CT_COMPILE_EXE "Compile the src/main.cpp executable" OFF)
option(AFF3CT_COMPILE_TESTS "Compile the tests in tests/main.cpp" OFF)

if (AFF3CT_COMPILE_EXE)
    message(STATUS "AFF3CT_COMPILE_EXE -ON")
endif()
if (AFF3CT_COMPILE_TESTS)
    message(STATUS "AFF3CT_COMPILE_TESTS -ON")
endif()

project (aff3ct_nand_flash)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify bin path
if (AFF3CT_COMPILE_EXE)
    set (EXECUTABLE_OUTPUT_PATH bin/)
endif()
if (AFF3CT_COMPILE_TESTS)
    set (EXECUTABLE_OUTPUT_PATH tests/)
endif()

# Fetch googletest
if (AFF3CT_COMPILE_TESTS)
    include(FetchContent)
    FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest
    GIT_TAG release-1.11.0
    )
    FetchContent_GetProperties(googletest)
    if(NOT googletest_POPULATED)
        FetchContent_Populate(googletest)
        add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BUILD_DIR})
    endif()
endif()

# Create the executable from sources
set(source_dir "")
if (AFF3CT_COMPILE_EXE)
    set (source_dir "src") 
    add_executable(aff3ct_nand_flash ${CMAKE_CURRENT_SOURCE_DIR}/${source_dir}/main.cpp)
endif()
if (AFF3CT_COMPILE_TESTS)
    set (source_dir "tests")
endif()

file(GLOB_RECURSE source_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
add_library(aff3ct_extension ${source_files}) #
target_include_directories(aff3ct_extension PUBLIC include/)


# Link with the "Threads library (required to link with AFF3CT after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# Link with google tests
if (AFF3CT_COMPILE_TESTS)
    file(GLOB_RECURSE test_source_files ${CMAKE_CURRENT_SOURCE_DIR}/tests/*)
    add_executable(aff3ct_nand_flash ${test_source_files})
    #add_library(aff3ct_tests ${test_source_files}) 
    #target_include_directories(aff3ct_tests PUBLIC tests/)
    #target_link_libraries(aff3ct_tests PRIVATE aff3ct::aff3ct-static-lib aff3ct_extension gtest_main)
    target_link_libraries(aff3ct_nand_flash PRIVATE gtest_main)
endif()
    

# Link with AFF3CT
set (AFF3CT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
find_package(AFF3CT CONFIG 2.3.2 REQUIRED)
target_link_libraries(aff3ct_extension PRIVATE aff3ct::aff3ct-static-lib)
target_link_libraries(aff3ct_nand_flash PRIVATE aff3ct::aff3ct-static-lib aff3ct_extension)