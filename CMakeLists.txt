cmake_minimum_required(VERSION 3.0.0)
project(TestProject VERSION 0.1.0)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#Add library
add_subdirectory(TestLibrary)
list(APPEND EXTRA_LIBS TestLibrary)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/TestLibrary")

#Build multiple executables
add_executable(test1 testcommand1.cpp)
add_executable(test2 testcommand2.cpp)

#Find directory
target_include_directories(test1 PUBLIC
"${PROJECT_BINARY_DIR}"
${EXTRA_INCLUDES}
)

target_include_directories(test2 PUBLIC
"${PROJECT_BINARY_DIR}"
${EXTRA_INCLUDES}
)

#Link library
target_link_libraries(test1 PUBLIC ${EXTRA_LIBS})
target_link_libraries(test2 PUBLIC ${EXTRA_LIBS})

#Testing
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
FetchContent_MakeAvailable(googletest)


enable_testing()

add_executable(
  testunit
  testunit.cpp
)
target_link_libraries(
  testunit
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(testunit)





