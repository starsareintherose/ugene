set(UGENE_PLUGIN_NAME GUITestBase)

include_directories(../../libs_3rdparty/QSpec/src)
include_directories(../../libs_3rdparty/QSpec/src/primitives)
include_directories(../../libs_3rdparty/QSpec/src/utils)
include_directories(src)

include(../../Plugin.cmake)

if(WIN32)
    set(CMAKE_CXX_FLAGS "-bigobj")
endif()

add_definitions(-DUNICODE)

if(UGENE_CUDA_DETECTED)
    add_definitions(-DSW2_BUILD_WITH_CUDA)
endif()