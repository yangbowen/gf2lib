cmake_minimum_required (VERSION 3.26)

project ("gf2lib" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS FALSE)

add_library("gf2lib" STATIC)
target_sources(
	"gf2lib"
	PRIVATE
	"clmul.h"
	"crc.h"
	"crc.cpp"
)

add_executable("gf2lib-test")
target_link_libraries("gf2lib-test" "gf2lib")
target_sources(
	"gf2lib-test"
	PRIVATE
	"clmul-test.cpp"
)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	target_compile_options(
		${PROJECT_NAME}
		PRIVATE
		"/W4"
	)
elseif(
	(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
)
	target_compile_options(
		${PROJECT_NAME}
		PRIVATE
		"-Wall"
		"-Wextra"
		"-Wno-c++98-compat"
		"-Wno-c++98-compat-pedantic"
		"-Wno-c++0x-compat"
		"-Wno-c++11-compat"
		"-Wno-c++14-compat"
		"-Wno-c++17-compat"
		"-Wno-c++20-compat"
	)
endif()
