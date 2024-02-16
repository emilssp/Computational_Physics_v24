cmake_minimum_required(VERSION 3.13)

project(quantum_mech_prg
        VERSION 0.1.0 
        LANGUAGES C CXX)

        set(prog_name
        mainEXE
        )
    
    set(lib_name
        QM_lib
        )
    
    add_executable(${prog_name} main.cpp)
    
    target_include_directories(
        ${prog_name}
        PRIVATE 
        ${PROJECT_SOURCE_DIR}/../src
        )
    
    if(TARGET ${lib_name})
       message(STATUS "it found target " ${lib_name})
       message(STATUS "Project source dir is " ${PROJECT_SOURCE_DIR})
    else()
       message(SEND_ERROR "Attempt to link to non-existent library: " ${lib_name})
    endif()
    
    target_link_libraries(
        ${prog_name}
        #${PROJECT_SOURCE_DIR}/../src/${lib_name}
        ${lib_name}
        )
    