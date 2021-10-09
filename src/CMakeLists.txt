set(HEADER_FILES
        ../include/XSecAna/_Hist.h
        ../include/XSecAna/Array.h
        ../include/XSecAna/Hist.h
)

set(SOURCES
        ./Hist.cpp
        ./_Hist.cpp
        ./Array.cpp
)


add_library(XSecAna SHARED
        ${HEADER_FILES}
        ${SOURCES}
)

target_link_libraries(XSecAna PRIVATE
        ${ROOT_LIBRARIES}
)
