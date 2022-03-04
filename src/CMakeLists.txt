set(HEADER_FILES
        ../include/XSecAna/SimpleEfficiency.h
        ../include/XSecAna/SimpleFlux.h
        ../include/XSecAna/SimpleSignalEstimator.h
        ../include/XSecAna/IUnfold.h
        ../include/XSecAna/CrossSection.h
        ../include/XSecAna/Systematic.h
        ../include/XSecAna/SimpleQuadSum.h
        ../include/XSecAna/SimpleQuadSumAsymm.h
        ../include/XSecAna/Fit/IFitter.h
        ../include/XSecAna/Fit/TemplateFitCalculator.h
        ../include/XSecAna/Fit/Minuit2TemplateFitter.h
        ../include/XSecAna/Fit/TemplateFitComponent.h
        ../include/XSecAna/TemplateFitSignalEstimator.h
        ../include/XSecAna/Fit/JointTemplateFitComponent.h
        ../include/XSecAna/JointTemplateFitSignalEstimator.h
)

set(SOURCES
        ./SimpleEfficiency.cpp
        ./SimpleFlux.cpp
        ./SimpleSignalEstimator.cpp
        ./CrossSection.cpp
        ./Systematic.cpp
        ./Fit/TemplateFitCalculator.cpp
        ./Fit/Minuit2TemplateFitter.cpp
        ./Fit/TemplateFitComponent.cpp
        ./TemplateFitSignalEstimator.cpp
        ./Fit/JointTemplateFitComponent.cpp
        ./JointTemplateFitSignalEstimator.cpp
)


add_library(XSecAna SHARED
        ${HEADER_FILES}
        ${SOURCES}
)

target_link_libraries(XSecAna PRIVATE
        ${ROOT_LIBRARIES}
        ROOT::Minuit2
)
