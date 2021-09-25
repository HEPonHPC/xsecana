//
// Created by Derek Doyle on 9/23/21.
//

#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/Fit/Minuit2TemplateFitter.h"

#include <Eigen/Dense>
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv) {
    auto nparams = 5;
    auto mask_param = 3;

    Eigen::VectorXd user_params(nparams);
    user_params << 1, 2, 3, 4, 5;
    Eigen::VectorXd minimizer_params(nparams);
    minimizer_params << 1, 2, 3, 4, 5;
    Eigen::VectorXd reduced_user_params(nparams);
    Eigen::VectorXd reduced_minimizer_params(nparams - 1);

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


    std::vector<int> dims = {4, 10};
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> signal_templates(dims[0], dims[1]);
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> background1_templates(dims[0], dims[1]);
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> background2_templates(dims[0], dims[1]);

    auto x = Eigen::Array<double, 1, -1>::LinSpaced(dims[1], 0, dims[1]);
    for (auto i = 0u; i < signal_templates.rows(); i++) {
        signal_templates.row(i) = (i + 1) * x + 0.5;
        background1_templates.row(i) = (i + 1) * x.reverse() + 0.5;
        background2_templates.row(i) = (i + 1) * Eigen::Matrix<double, 1, -1>::Ones(dims[1]);
    }

    std::vector<Eigen::Array<double, 1, -1>> templates{
            signal_templates.reshaped<Eigen::RowMajor>().transpose(),
            background1_templates.reshaped<Eigen::RowMajor>().transpose(),
            background2_templates.reshaped<Eigen::RowMajor>().transpose(),
    };
    Eigen::VectorXd total = templates[0] + templates[1] + templates[2];

    Eigen::MatrixXd inverse_covariance = Eigen::MatrixXd::Identity(dims[0] * dims[1],
                                                                   dims[0] * dims[1]);

    inverse_covariance *= (1 / total.array()).matrix().asDiagonal();

    auto fit_calc = new fit::TemplateFitCalculator(templates, dims, inverse_covariance);

    user_params = Eigen::RowVectorXd::Ones(templates.size() * dims[0]);
    minimizer_params = Eigen::RowVectorXd::Ones(templates.size() * dims[0] - 1);

    assert((total - fit_calc->Predict(user_params)).isZero(0));

    user_params(user_params.size() - 1) = 4;
    fit_calc->FixTemplate(fit_calc->GetNTemplates() - 1, 4);
    assert((fit_calc->U(minimizer_params) - fit_calc->Predict(user_params)).isZero());

    // make sure we can fix the same template to a different value
    user_params(user_params.size() - 1) = 5;
    fit_calc->FixTemplate(fit_calc->GetNTemplates() - 1, 5);
    assert((fit_calc->U(minimizer_params) - fit_calc->Predict(user_params)).isZero());

    // now release a template
    fit_calc->ReleaseTemplate(fit_calc->GetNTemplates() - 1);
    user_params = Eigen::RowVectorXd::Ones(templates.size() * dims[0]);
    assert((total - fit_calc->U(user_params)).isZero(0));


    user_params(4) = 1.5;

    fit::Minuit2TemplateFitter fitter(fit_calc, 3);
    fitter.SetPrintLevel(1);
    auto result = fitter.Fit(fit_calc->Predict(user_params));

    std::cout << result.params << std::endl;
    std::cout << result.fun_val << std::endl;
    std::cout << result.fun_calls << std::endl;
}


