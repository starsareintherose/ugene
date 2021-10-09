# add cuda dependency
if (UGENE_CUDA_DETECTED)
    set(UGENE_PLUGIN_NAME cuda_support)

    include(../../Plugin.cmake)

    include_directories(${UGENE_CUDA_INC_DIR})
endif()