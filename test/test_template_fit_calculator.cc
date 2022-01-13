//
// Created by Derek Doyle on 9/23/21.
//


#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/Fit/Minuit2TemplateFitter.h"
#
#include <Eigen/Dense>
#include <iostream>

#include "TFile.h"

using namespace xsec;
using namespace xsec::fit;
int main(int argc, char ** argv) {
    auto nparams = 5;
    auto mask_param = 3;

    Vector user_params(nparams);
    user_params << 1, 2, 3, 4, 5;
    Vector minimizer_params(nparams);
    minimizer_params << 1, 2, 3, 4, 5;
    Vector reduced_user_params(nparams);
    Vector reduced_minimizer_params(nparams - 1);

    auto reduced_idx = 0;
    for (auto i = 0u; i < nparams; i++) {
        if (i == mask_param) {
            reduced_user_params(i) = 0;
        } else {
            reduced_user_params(i) = user_params(i);
            reduced_minimizer_params(reduced_idx) = minimizer_params(i);
            reduced_idx++;
        }
    }

    fit::detail::ParamMap param_map(nparams);
    for (auto i = 0u; i < nparams; i++) {
        assert(!param_map.IsParamMasked(i));
    }
    assert((user_params - param_map.ToUserParams(minimizer_params)).isZero(0));
    assert((minimizer_params - param_map.ToMinimizerParams(user_params)).isZero(0));

    param_map.MaskTemplate(mask_param);
    assert(param_map.IsParamMasked(mask_param));
    for (auto i = 0u; i < nparams; i++) {
        if (i == mask_param) assert(param_map.IsParamMasked(i));
        else
            assert(!param_map.IsParamMasked(i));
    }
    assert((reduced_user_params - param_map.ToUserParams(reduced_minimizer_params)).isZero(0));
    assert((reduced_minimizer_params - param_map.ToMinimizerParams(user_params)).isZero(0));

    param_map.UnmaskTemplate(mask_param);
    for (auto i = 0u; i < nparams; i++) {
        assert(!param_map.IsParamMasked(i));
    }
    assert((user_params - param_map.ToUserParams(minimizer_params)).isZero(0));
    assert((minimizer_params - param_map.ToMinimizerParams(user_params)).isZero(0));

    param_map.MaskTemplate(mask_param);
    for (auto i = 0u; i < nparams; i++) {
        if (i == mask_param) assert(param_map.IsParamMasked(i));
        else
            assert(!param_map.IsParamMasked(i));
    }
    assert((reduced_user_params - param_map.ToUserParams(reduced_minimizer_params)).isZero(0));
    assert((reduced_minimizer_params - param_map.ToMinimizerParams(user_params)).isZero(0));

    std::vector<double> std_vector;
    for(auto i = 0u; i < 10; i++) std_vector.push_back(i);
    Vector from_std = fit::detail::STDToEigen(std_vector);
    for(auto i = 0u; i < std_vector.size(); i++) {
        assert(from_std[i] == std_vector[i]);
    }

    std::vector<int> dims = {4, 10};
    Matrix signal_templates(dims[0], dims[1]);
    Matrix background1_templates(dims[0], dims[1]);
    Matrix background2_templates(dims[0], dims[1]);

    auto x = Array::LinSpaced(dims[1], 0, dims[1]);
    for (auto i = 0u; i < signal_templates.rows(); i++) {
        signal_templates.row(i) = (i + 1) * x + 0.5;
        background1_templates.row(i) = (i + 1) * x.reverse() + 0.5;
        background2_templates.row(i) = (i + 1) * Vector::Ones(dims[1]);
    }

    std::map<std::string, IReducedTemplateComponent*> templates = {
            {"a", new ReducedTemplateComponent(new ReducedComponent(signal_templates.reshaped().transpose(), dims[0], dims[1]))},
            {"b", new ReducedTemplateComponent(new ReducedComponent(background1_templates.reshaped().transpose(), dims[0], dims[1]))},
            {"c", new ReducedTemplateComponent(new ReducedComponent(background2_templates.reshaped().transpose(), dims[0], dims[1]))},
    };
    Vector total = templates["a"]->GetNominal()->GetArray() + templates["b"]->GetNominal()->GetArray() +
                   templates["c"]->GetNominal()->GetArray();

    Matrix inverse_covariance = Matrix::Identity(dims[0] * dims[1],
                                                 dims[0] * dims[1]);

    inverse_covariance *= (1 / total.array()).matrix().asDiagonal();

    auto fit_calc = new TemplateFitCalculator(ReducedComponentCollection(templates), dims, inverse_covariance);

    user_params = Vector::Ones(templates.size() * dims[0]);

    assert((total - fit_calc->Predict(user_params)).isZero(0));

    // now release a template
    fit_calc->ReleaseTemplate(fit_calc->GetNTemplates() - 1);
    user_params = Vector::Ones(templates.size() * dims[0]);
    assert((total - fit_calc->Predict(user_params)).isZero(0));

    user_params(4) = 1.5;

    fit::Minuit2TemplateFitter fitter(3);

    fitter.SetPrintLevel(0);
    auto result = fitter.Fit(fit_calc, fit_calc->Predict(user_params));

    assert(result.fun_val == fit_calc->Chi2(fit_calc->ToUserParams(result.params),
                                            fit_calc->Predict(user_params)));
}


