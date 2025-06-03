# tools/CMakeLists.txt
add_library(ecoli_parser STATIC
    ${PROJECT_SOURCE_DIR}/src/ecoli_parser.cpp)
target_include_directories(ecoli_parser PUBLIC
    ${PROJECT_SOURCE_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(CURL REQUIRED libcurl)
find_package(ZLIB REQUIRED)

target_link_libraries(ecoli_parser PRIVATE ${CURL_LIBRARIES} ${ZLIB_LIBRARIES})
target_compile_options(ecoli_parser PRIVATE ${CURL_CFLAGS_OTHER})

add_executable(bamboofilter_ecoli
    ${PROJECT_SOURCE_DIR}/src/run_ecoli.cpp)
target_link_libraries(bamboofilter_ecoli
    ecoli_parser
    header  
    hash    
    pthread
    ${CURL_LIBRARIES} ${ZLIB_LIBRARIES})
target_include_directories(bamboofilter_ecoli PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src 
    ${CURL_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
