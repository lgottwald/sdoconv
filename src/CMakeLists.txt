cmake_minimum_required(VERSION 2.8)

find_package(libsdo REQUIRED)

ADD_EXECUTABLE(sdoconv
	GamsGenerator.cpp
	Main.cpp
	Escape.cpp
	)

include_directories(${libsdo_INCLUDE_DIRS})
set_target_properties(sdoconv PROPERTIES COMPILE_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Wno-unused-parameter")

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

TARGET_LINK_LIBRARIES(sdoconv ${libsdo_LIBRARIES} ${Boost_LIBRARIES} )

INSTALL(TARGETS sdoconv RUNTIME DESTINATION bin)
