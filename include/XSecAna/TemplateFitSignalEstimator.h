#pragma once
#include "XSecAna/Systematic.h"
#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/Fit/IFitter.h"

#include "TCanvas.h"

namespace xsec {
    struct TemplateFitResult {
        double fun_val;
        TH1 * signal_params;
        TH1 * signal_params_error_up;
        TH1 * signal_params_error_down;
        std::map<std::string, TH1*> background_params;
        std::map<std::string, TH1*> background_params_error_up;
        std::map<std::string, TH1*> background_params_error_down;
        TH2D * covariance;
        unsigned int fun_calls;
    };

    class TemplateFitSignalEstimator : public IEigenSignalEstimator {
    public:

        TemplateFitSignalEstimator(const TH1 * signal_template,
                                   const std::map<std::string, const TH1*> & background_templates,
                                   const std::map<std::string, Systematic<TH1>> & systematics,
                                   const TH1 * mask = 0,
                                   fit::TemplateFitCalculator * fit_calc = 0);

        void SetFitter(fit::IFitter * fitter);
        fit::IFitter * GetFitter() const;
        fit::IFitCalculator * GetFitCalc() const;

        TemplateFitResult Fit(const TH1 * data, int nrandom_seeds=-1) const;
        TemplateFitResult Fit(const TH1 * data, fit::IFitter * fitter, int nrandom_seeds=-1);

        virtual TH1 * Background(const TH1 * data) const override;

        virtual TH1 * Signal(const TH1 * data) const override;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const override;

        void FixComponent(const std::string & template_name, const double & val=1);
        void ReleaseComponent(const std::string & template_name);

        /// \brief given params, predict the number of signal events
        /// in each template bin by weighing templates
        TH1 * PredictTotal(const TH1* signal_params,
                           const std::map<std::string, TH1*> & bkgd_params) const;
        TH1 * PredictSignal(const TH1 * signal_params) const;
        TH1 * PredictBackground(const std::string & background_label,
                                const TH1 * bkgd_params) const;

        /// \brief given params, predict the number of signal events
        /// in each analysis bin by weighing templates and integrated template bins
        TH1 * PredictProjectedTotal(const TH1* signal_params,
                                    const std::map<std::string, TH1*> & bkgd_params) const;
        TH1 * PredictProjectedSignal(const TH1 * signal_params) const;
        TH1 * PredictProjectedBackground(const std::string & background_label,
                                         const TH1 * bkgd_params) const;
        double Chi2(const TH1 * data,
                    const TH1 * signal_params,
                    const std::map<std::string, TH1*> & bkgd_params) const;
        TH2D * GetTotalCovariance() const;
        TH2D * GetCovariance(const std::string & systematic_name) const;
        TH2D * GetInverseCovariance() const;
        TH1D * GetReducedSignalTemplate() const;
        TH1D * GetReducedBackgroundTemplate(const std::string & bkgd_label) const;
        TH1D * GetReducedTotalTemplate() const;

        TCanvas * DrawParameterCovariance(const TemplateFitResult & fit_result) const;
        TCanvas * DrawParameterCorrelation(const TemplateFitResult & fit_result) const;

        std::vector<std::string> GetBackgroundLabels() const;
        std::vector<std::string> GetSystematicLabels() const;

    private:
        void _draw_covariance_helper(TCanvas * c, TH1 * mat, const TemplateFitResult & fit_result) const;
        TemplateFitResult _template_fit_result(const fit::FitResult & fit_result) const;
        TemplateFitResult _fit(const TH1 * data, const std::vector<Vector> & seeds) const;
        TH1 * _mask_and_flatten(const TH1 * mask, const TH1 * templ) const;
        TH1 * _mask_and_flatten(const TH1 * templ) const;
        TH1 * _to_template_binning(const Array & reduced_templates) const;
        Array _predict_component(const TH1 * component_templates, const TH1 * params) const;
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;

        fit::Vector ToCalculatorParams(const TH1 * signal_params,
                                       const std::map<std::string, TH1*> & bkgd_params) const;

        void ToUserParams(const fit::Vector & calc_params,
                          TH1 * signal_params,
                          std::map<std::string, TH1*> & bkgd_params) const;

        TH1 * fSignalTemplate;
        std::map<std::string, TH1*> fBackgroundTemplates;
        std::map<std::string, Systematic<TH1>> fSystematics;
        TH1 * fTotalTemplate;

        std::vector<int> fDims;

        std::map<std::string, bool> fIsFreeTemplate;

        fit::TemplateFitCalculator * fFitCalc;
        std::map<std::string, Matrix> fCovarianceMatrices;
        Matrix fTotalCovariance;
        Matrix fInverseCovariance;

        std::map<std::string, int> fComponentLabelIdxMap;

        const TH1 * fMask;

        root::TH1Props fProjectPredictionProps;
        root::TH1Props fPredictionProps;
        fit::detail::ParamMap fOuterBinMap;
        fit::detail::ParamMap fInnerBinMap;

        fit::IFitter * fFitter = 0;
    };
}