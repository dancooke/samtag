set(CXX_OPTIMIZATION_FLAGS -march=${COMPILER_ARCHITECTURE})

add_executable(samtag main.cpp)
target_compile_features(samtag PRIVATE cxx_std_20)
target_compile_options(samtag PRIVATE ${CXX_OPTIMIZATION_FLAGS})

find_package (HTSlib 1.14 REQUIRED)
target_include_directories (samtag PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries (samtag ${HTSlib_LIBRARIES})

check_ipo_supported(RESULT ipo_supported OUTPUT output)
if(ipo_supported)
    message("-- IPO is supported!")
    set_property(TARGET samtag PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()

install(TARGETS samtag DESTINATION ${CMAKE_INSTALL_PREFIX})
