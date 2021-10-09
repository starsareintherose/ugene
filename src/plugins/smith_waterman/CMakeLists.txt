set(UGENE_PLUGIN_NAME smith_waterman)

#add_definitions(-DSW2_BUILD_WITH_SSE2)


if(UGENE_CUDA_DETECTED)
    find_package(CUDAToolkit)

    #LIBS += -L$$UGENE_CUDA_LIB_DIR -lcudart
    set(UGENE_ADDITIONAL_LIBS CUDA::cudart)
    include_directories(${UGENE_CUDA_INC_DIR})

    set(SW2_CUDA_FILES src/sw_cuda.cu)
endif()

include(../../Plugin.cmake)