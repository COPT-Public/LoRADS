cmake_minimum_required(VERSION 3.21)
project(LoRADS_v_2_0_1-alpha C)
set(CMAKE_C_STANDARD 11)

# set(MODE "STATIC") # set(MODE "DYNAMIC")
set(MODE "DYNAMIC")
if(MODE STREQUAL "STATIC")
    if (APPLE)
        set(CMAKE_OSX_ARCHITECTURES x86_64)
        set(MKL_HOME /opt/intel/oneapi/mkl/2023.2.2/lib)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(ARPACK REQUIRED arpack)
    endif()


    if (UNIX AND NOT APPLE)
        set(MKL_HOME /usr/lib/x86_64-linux-gnu)
    endif()
   set(CMAKE_BUILD_TYPE Release)
    # set(CMAKE_BUILD_TYPE Debug)
    add_compile_definitions(MAC_INT64)
    include_directories(externals
            src_semi
            src_semi/data
            src_semi/lorads_alg
            src_semi/linalg
            src_semi/io
    )
    #    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    file(GLOB SOURCES src_semi/*.c src_semi/data/*.c src_semi/linalg/*.c src_semi/lorads_alg/*.c externals/*.c src_semi/io/*.c)
    file(GLOB HEADERS src_semi/*.h src_semi/data/*.h src_semi/linalg/*.h src_semi/lorads_alg/*.h externals/*.h src_semi/io/*.h)
    add_executable(LoRADS_v_2_0_1-alpha
            ${SOURCES}
            ${HEADERS}
    )
    if (APPLE)
        target_link_libraries(LoRADS_v_2_0_1-alpha
                ${MKL_HOME}/libmkl_intel_lp64.a
                ${MKL_HOME}/libmkl_sequential.a
                ${MKL_HOME}/libmkl_core.a
                m
                ${ARPACK_LIBRARIES}
                ${ARPACK_LDFLAGS}
        )
    endif()

    if (UNIX AND NOT APPLE)
        target_link_libraries(LoRADS_v_2_0_1-alpha
                -Wl,--start-group
                ${MKL_HOME}/libmkl_intel_lp64.a
                ${MKL_HOME}/libmkl_sequential.a
                ${MKL_HOME}/libmkl_core.a
                -Wl,--end-group
                pthread
                m
                dl
        )
    endif()
else()
    # dynamic
    if (APPLE)
        find_library(ACCELERATE_FRAMEWORK Accelerate)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(ARPACK REQUIRED arpack)
        add_compile_definitions(MAC_INT64)
    endif()


    if (UNIX AND NOT APPLE)
        find_library(LAPACK_LIBRARIES NAMES lapack PATHS /usr/lib /usr/local/lib  /lib/x86_64-linux-gnu/)
        find_library(BLAS_LIBRARIES NAMES blas PATHS /usr/lib /usr/local/lib  /lib/x86_64-linux-gnu/)
        add_compile_definitions(INT32)
        if(NOT LAPACK_LIBRARIES)
            message(FATAL_ERROR "LAPACK library not found")
        endif()

        if(NOT BLAS_LIBRARIES)
            message(FATAL_ERROR "BLAS library not found")
        endif()
        add_compile_definitions(UNDER_BLAS)
    endif()
   set(CMAKE_BUILD_TYPE Release)
    #  set(CMAKE_BUILD_TYPE Debug)
    #    add_compile_definitions(FIX_INI_POINT)
    # add_compile_definitions(DUAL_U_V)
    # add_compile_definitions(INT32)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    include_directories(externals
            src_semi
            src_semi/data
            src_semi/lorads_alg
            src_semi/linalg
            src_semi/io
    )
    file(GLOB SOURCES src_semi/*.c src_semi/data/*.c src_semi/linalg/*.c src_semi/lorads_alg/*.c externals/*.c src_semi/io/*.c)
    file(GLOB HEADERS src_semi/*.h src_semi/data/*.h src_semi/linalg/*.h src_semi/lorads_alg/*.h externals/*.h src_semi/io/*.h)
    add_executable(LoRADS_v_2_0_1-alpha
            ${SOURCES}
            ${HEADERS}
    )
    if (APPLE)
        target_link_libraries(LoRADS_v_2_0_1-alpha
                ${ACCELERATE_FRAMEWORK}
                ${ARPACK_LIBRARIES}
                ${ARPACK_LDFLAGS}
        )
    endif()

    if (UNIX AND NOT APPLE)
        target_link_libraries(LoRADS_v_2_0_1-alpha
                PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                m
                arpack)
    endif()
endif()
