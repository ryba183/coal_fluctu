cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(coal_fluctu CXX C)


#search for libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()


add_executable(coal_fluctu coal_fluctu.cpp)

target_compile_features(coal_fluctu PUBLIC cxx_std_11)

target_link_libraries(coal_fluctu
  PRIVATE
    clphxx::cloudphxx_lgrngn
)

