#pragma once

#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/Hist.h"
#include "XSecAna/ISignalEstimator.h"

namespace xsec {
    template<class HistType = HistXd>
    class TemplateFitSignalEstimator : public ISignalEstimator<HistType> {
    public:
        TemplateFitSignalEstimator(HistType signal_template,
                                   std::map<std::string, HistType> background_templates,
                                   std::vector<int> dims,
                                   fit::TemplateFitCalculator<HistType::scalar_type, HistType::cols_type> * fit_calc)
                : fTemplates(templates), fDims(dims), fFitCalc(fit_calc)
        {
            for(auto temp_it : fTemplates) {
                fIsFreeTemplate[temp_it->first] = true;
            }
        }

        virtual HistXd Eval(const HistType & data);

        virtual const HistXd & Background(const HistType & data);

        virtual const HistXd & Signal(const HistType & data) ;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const;

        void FixTemplate(std::string template_name)
        { fTemplates.at(template_name).second = true;}

        /// \brief given params, predict the number of signal events
        /// by weighing templates and integrating over the last dimension
        HistXd & Predict(const std::vector<double> & params, const HistType & data) const;

    private:

        std::vector<int> GetSignalTemplateIdxs();


        HistType fSignalTemplate;
        std::map<std::string, HistType> fBackgroundTemplates;
        std::vector<int> fDims;
        std::map<std::string, bool> fIsFreeTemplate;

        TemplateFitCalculator<HistType> * fFitCalc;

    };

    template<class HistType>
    HistType
    TemplateFitSignalEstimator<HistType>::
    Eval(const HistType & data) {
        auto [params, errors] = fFitCalc->Fit(this, data);
        auto prediction = this->Predict(params, data);




    }
}