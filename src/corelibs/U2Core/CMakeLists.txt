project(U2Core)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED Core Gui Widgets Script Network Xml)

add_definitions(
        -DBUILDING_U2CORE_DLL
)

include_directories(src)
include_directories(../../include)
include_directories(../../libs_3rdparty/sqlite3)
include_directories(../../libs_3rdparty/sqlite3/src)

file(GLOB_RECURSE SRCS src/*.cpp)
file(GLOB_RECURSE HDRS src/*.h)

add_library(U2Core SHARED ${HDRS} ${SRCS})

if(APPLE)
    target_link_libraries(U2Core "-framework CoreFoundation")
endif()

target_link_libraries(U2Core
        Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Network Qt5::Script Qt5::Xml
        zlib ugenedb)
