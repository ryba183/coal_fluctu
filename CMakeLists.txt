cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(coal_fluctu CXX C)


#search for libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

set(CMAKE_CXX_FLAGS ${libcloudphxx_CXX_FLAGS_RELEASE})

add_executable(coal_fluctu coal_fluctu.cpp)

target_link_libraries(coal_fluctu ${libcloudphxx_LIBRARIES})
target_include_directories(coal_fluctu PUBLIC ${libcloudphxx_INCLUDE_DIRS})

