#pragma once
#include "XSecAna/Systematic.h"
#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/Fit/IFitter.h"

namespace xsec {
    class TemplateFitSignalEstimator : public IEigenSignalEstimator {
    public:
        TemplateFitSignalEstimator(const TH1 * signal_template,
                                   const std::map<std::string, const TH1*> & background_templates,
                                   const std::map<std::string, Systematic<TH1>> & systematics,
                                   const TH1 * mask = 0,
                                   fit::TemplateFitCalculator * fit_calc = 0);
        virtual TH1 * Background(const TH1 * data) const override;

        virtual TH1 * Signal(const TH1 * data) const override;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const override;

        void FixTemplate(const std::string & template_name);

        /// \brief given params, predict the number of signal events
        /// in each template bin by weighing templates
        TH1 * Predict(const TH1* signal_params,
                      const std::map<std::string, TH1*> & bkgd_params) const;

        /// \brief given params, predict the number of signal events
        /// in each analysis bin by weighing templates and integrated template bins
        TH1 * PredictProjected(const TH1* signal_params,
                               const std::map<std::string, TH1*> & bkgd_params) const;

        //std::pair<TH1*, std::map<std::string, TH1*>>
        //PredictComponents(const Matrix & signal_params,
        //                  const std::map<std::string, Matrix> & bkgd_params) const;
        void SetMask(const Array & mask);

        TH2D * GetTotalCovariance() const;
        TH2D * GetCovariance(const std::string & systematic_name) const;
        TH2D * GetInverseCovariance() const;
        TH2D * GetReducedSignalTemplate() const;
        TH2D * GetReducedBackgroundTemplate(const std::string & bkgd_label) const;
        TH2D * GetReducedTotalTemplate() const;
        const Systematic<TH1> & GetReducedSystematic(const std::string & systematic_label) const;


    private:
        TH1 * _mask_and_flatten(const TH1 * mask, const TH1 * templ) const;
        const TH1 * _transpose(const TH1 * h) const;
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;

        fit::Vector ToCalculatorParams(const TH1 * signal_params,
                                       const std::map<std::string, TH1*> & bkgd_params) const;

        void ToUserParams(const fit::Vector & calc_params,
                          TH1 * signal_params,
                          std::map<std::string, TH1*> & bkgd_params) const;

        std::vector<int> GetSignalTemplateIdxs();


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

        std::map<int, std::string> fComponentLabelIdxMap;


        const TH1 * fMask;

        root::TH1Props fProjectPredictionProps;
        fit::detail::ParamMap fParamMap;
    };
}