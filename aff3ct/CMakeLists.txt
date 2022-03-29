cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0054 NEW)

project (aff3ct_nand_flash)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify bin path
set (EXECUTABLE_OUTPUT_PATH bin/)

# Create the executable from sources
add_executable(aff3ct_nand_flash ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

# Create library for extension code
# Find all subdirectories macro
#macro(subdirlist result curdir)
#  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
#  set(dirlist "")
#  foreach(child ${children})
#    if(IS_DIRECTORY ${curdir}/${child})
#      list(APPEND dirlist ${child})
#    endif()
#  endforeach()()
#  set(${result} ${dirlist})
#endmacro()

#subdirlist(INCLUDE_SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
#subdirlist(SRC_SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB_RECURSE source_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
add_library(aff3ct_extension ${source_files}) #
target_include_directories(aff3ct_extension PUBLIC include/)


# Link with the "Threads library (required to link with AFF3CT after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Link with AFF3CT
set (AFF3CT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
find_package(AFF3CT CONFIG 2.3.2 REQUIRED)
target_link_libraries(aff3ct_extension PRIVATE aff3ct::aff3ct-static-lib)
target_link_libraries(aff3ct_nand_flash PRIVATE aff3ct::aff3ct-static-lib aff3ct_extension)