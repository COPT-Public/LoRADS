cmake_minimum_required(VERSION 3.21)
project(LoRADS_v_1_0_0-alpha C)
set(CMAKE_C_STANDARD 11)

set(MODE "STATIC") # set(MODE "DYNAMIC")

if(MODE STREQUAL "STATIC")
    if (APPLE)
        set(CMAKE_OSX_ARCHITECTURES x86_64)
        set(MKL_HOME /opt/intel/oneapi/mkl/2023.1.0/lib)
    endif()


    if (UNIX AND NOT APPLE)
        set(MKL_HOME /usr/lib/x86_64-linux-gnu)
    endif()
    set(CMAKE_BUILD_TYPE Release)
    include_directories(externals
                        interface
                        src
                        )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    file(GLOB SOURCES src/*.c externals/*.c)
    file(GLOB HEADERS src/*.h externals/*.h)
    add_executable(LoRADS_v_1_0_0-alpha
                ${SOURCES}
                ${HEADERS}
                )
    if (APPLE)
        target_link_libraries(LoRADS_v_1_0_0-alpha
            ${MKL_HOME}/libmkl_intel_lp64.a
            ${MKL_HOME}/libmkl_sequential.a
            ${MKL_HOME}/libmkl_core.a
            m
        )
    endif()

    # if (UNIX AND NOT APPLE)
    # target_link_libraries(LoRADS_v_1_0_0-alpha
    #     -Wl,--start-group
    #     ${MKL_HOME}/libmkl_intel_lp64.a
    #     ${MKL_HOME}/libmkl_intel_thread.a
    #     ${MKL_HOME}/libmkl_core.a
    #     -Wl,--end-group
    #     iomp5
    #     pthread
    #     m
    #     dl
    # )
    # endif()

    if (UNIX AND NOT APPLE)
    target_link_libraries(LoRADS_v_1_0_0-alpha
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
    endif()


    if (UNIX AND NOT APPLE)
        find_package(LAPACK REQUIRED)
    endif()
    set(CMAKE_BUILD_TYPE Release)
    include_directories(externals
                    interface
                    src
                    )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    file(GLOB SOURCES src/*.c externals/*.c)
    file(GLOB HEADERS src/*.h externals/*.h)
    add_executable(LoRADS_v_1_0_0-alpha
            ${SOURCES}
            ${HEADERS}
            )
    if (APPLE)
    target_link_libraries(LoRADS_v_1_0_0-alpha
                        ${ACCELERATE_FRAMEWORK}
                        )
    endif()

    if (UNIX AND NOT APPLE)
    target_link_libraries(LoRADS_v_1_0_0-alpha
                            LAPACK::LAPACK
                            m
    )
    endif()
endif()