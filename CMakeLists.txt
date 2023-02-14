cmake_minimum_required(VERSION 3.9.0)

project(LACTRec VERSION 1.0)
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/compiled")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_EXPORT_COMPILE_COMMANDS 1 CACHE BOOL "for clang" FORCE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "  -O2  ")

find_package(ROOT  CONFIG REQUIRED COMPONENTS Minuit)
include("${ROOT_USE_FILE}")
root_generate_dictionary(LactRecEvent ${PROJECT_SOURCE_DIR}/include/LACTRecEvent.h 
		LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)


add_executable(LACTRec)
target_sources(LACTRec PUBLIC ${PROJECT_SOURCE_DIR}/src/LACTRec.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Reconstruction.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp ${PROJECT_SOURCE_DIR}/src/LACTRec.cpp ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_include_directories(LACTRec PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(LACTRec PRIVATE  ${ROOT_LIBRARIES})

add_executable(LACTAna)
target_include_directories(LACTAna PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_sources(LACTAna PUBLIC ${PROJECT_SOURCE_DIR}/src/Ana_Direction.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp  ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_link_libraries(LACTAna PRIVATE ${ROOT_LIBRARIES})



add_executable(MakeLookup)
target_include_directories(MakeLookup PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_sources(MakeLookup PUBLIC ${PROJECT_SOURCE_DIR}/src/MakeLookup.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp  ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_link_libraries(MakeLookup PRIVATE ${ROOT_LIBRARIES})
add_executable(DrawShape)
target_include_directories(DrawShape PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_sources(DrawShape PUBLIC ${PROJECT_SOURCE_DIR}/src/Draw_Shape.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp  ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_link_libraries(DrawShape PRIVATE ${ROOT_LIBRARIES})
add_executable(DrawQ)
target_include_directories(DrawQ PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_sources(DrawQ PUBLIC ${PROJECT_SOURCE_DIR}/src/ComputeQ.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp  ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_link_libraries(DrawQ PRIVATE ${ROOT_LIBRARIES})
