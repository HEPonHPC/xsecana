#pragma once
#include "XSecAna/Systematic.h"
#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/Fit/IFitter.h"

#include "TCanvas.h"

namespace xsec {
    struct TemplateFitResult {
        double fun_val;
        std::map<std::string, TH1*> component_params;
        std::map<std::string, TH1*> component_params_error_up;
        std::map<std::string, TH1*> component_params_error_down;
        TH2D * covariance;
        unsigned int fun_calls;

        TemplateFitResult Clone() const;

        void SaveTo(TDirectory * dir, const std::string & name) const;
        static std::unique_ptr<TemplateFitResult> LoadFrom(TDirectory * dir, const std::string & subdir);
    };

    class TemplateFitSignalEstimator {//}; : public IEigenSignalEstimator {
    public:
        TemplateFitSignalEstimator(const fit::TemplateFitSample & sample,
                                   const TH1 * mask = 0);

        void SetFitter(fit::IFitter * fitter);
        fit::IFitter * GetFitter() const;
        fit::IFitCalculator * GetFitCalc() const;

        TemplateFitResult Fit(const std::shared_ptr<TH1> data, int nrandom_seeds=-1) const;
        TemplateFitResult Fit(const std::shared_ptr<TH1> data, fit::IFitter * fitter, int nrandom_seeds=-1);

        TemplateFitResult Fit(const std::shared_ptr<TH1> data, const std::map<std::string, TH1*> & seed) const;
        TemplateFitResult Fit(const std::shared_ptr<TH1> data, fit::IFitter * fitter, const std::map<std::string, TH1*> & seed);

        void SaveTo(TDirectory * dir, const std::string & name) const;

        void FixComponent(const std::string & template_name, const double & val=1);
        void ReleaseComponent(const std::string & template_name);

        const fit::IUserTemplateComponent * GetUserComponent(const std::string & component_label) const
        { return fUserComponents.GetComponents().at(component_label); }

        /// \brief given params, predict the number of signal events
        /// in each template bin by weighing templates
        TH1 * PredictTotal(const std::map<std::string, TH1*> & params) const;
        TH1 * PredictComponent(const std::string & component_label, const TH1 * params) const;
        TH1 * NominalTotal() const;

        /// \brief given params, predict the number of signal events
        /// in each analysis bin by weighing templates and integrated template bins
        TH1 * PredictProjectedTotal(const std::map<std::string, TH1*> & params) const;
        TH1 * PredictProjectedComponent(const std::string & component_label, const TH1 * params) const;
        TH1 * NominalProjectedTotal() const;

        double Chi2(const std::shared_ptr<TH1> data,
                    const std::map<std::string, TH1*> & params = std::map<std::string, TH1*>()) const;
        TH2D * GetTotalCovariance() const;
        TH2D * GetCovariance(const std::string & systematic_name) const;
        TH2D * GetInverseCovariance() const;

        TCanvas * DrawParameterCovariance(const TemplateFitResult & fit_result) const;
        TCanvas * DrawParameterCorrelation(const TemplateFitResult & fit_result) const;

        std::vector<std::string> GetSystematicLabels() const;
        const fit::UserComponentCollection & GetUserComponentCollection() { return fUserComponents; }
        const fit::ReducedComponentCollection & GetReducedComponentCollection() { return fReducedComponents; }
        const std::map<std::string, const fit::IUserTemplateComponent*> & GetUserComponents() { return fUserComponents.GetComponents(); }
        const std::map<std::string, const fit::IReducedTemplateComponent*> & GetReducedComponents() { return fReducedComponents.GetComponents(); }

        const fit::IReducedTemplateComponent * GetReducedComponent(const std::string & component_label)
        { return fReducedComponents.GetComponents().at(component_label); }
        
        Systematic<TH1> PrefitComponentUncertainty(const std::string & component_label) const;
        TH1 * PostfitTotalUncertainty(const TemplateFitResult & fit_result) const;

        const TH1 * GetMask() const { return fReducer.GetMask(); }
        const fit::ComponentReducer & GetReducer() const { return fReducer; }

        fit::Vector ToCalculatorParams(const std::map<std::string, TH1*> & params) const;
        fit::Vector ToCalculatorParamsComponent(const TH1 * component_params) const;

        void ToUserParams(const fit::Vector & calc_params,
                          std::map<std::string, TH1 *> & params) const;
        TH1 * ToUserParamsComponent(const fit::Vector & calc_params) const;

    protected:
        bool _is_component_fixed(std::string label) const;
        void _draw_covariance_helper(TCanvas * c, TH1 * mat, const TemplateFitResult & fit_result) const;
        TemplateFitResult _template_fit_result(const fit::FitResult & fit_result) const;
        TemplateFitResult _fit(const std::shared_ptr<TH1> data, const std::vector<Vector> & seeds) const;
        TH1 * _to_template_binning(const Array & reduced_templates) const;
        TH1 * _project_prediction(const Array & prediction) const;



        fit::UserComponentCollection fUserComponents;
        fit::ReducedComponentCollection fReducedComponents;
        std::map<std::string, int> fComponentIdx;
        std::map<std::string, Systematic<TH1>> fSystematics;
        TH1 * fTotalTemplate = nullptr;

        std::map<std::string, const fit::IUserTemplateComponent*> fFixedUserComponents;
        fit::TemplateFitCalculator * fFitCalc;

        std::map<std::string, TH1*> fCovarianceMatrices;
        TH1 * fTotalCovariance;
        TH1 * fInvTotalCovariance;

        fit::ComponentReducer fReducer;

        root::TH1Props fProjectPredictionProps;
        root::TH1Props fPredictionProps;
        fit::detail::ParamMap fOuterBinMap;
        fit::detail::ParamMap fInnerBinMap;

        fit::IFitter * fFitter = 0;
    };
}