//
// Created by Derek Doyle on 9/23/21.
//

#include "XSecAna/Fit/TemplateFitCalculator.h"
#include "XSecAna/Fit/Minuit2TemplateFitter.h"

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

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

    std::vector<double> std_vector;
    for(auto i = 0u; i < 10; i++) std_vector.push_back(i);
    Eigen::VectorXd from_std = fit::detail::STDToEigen(std_vector);
    for(auto i = 0u; i < std_vector.size(); i++) {
        assert(from_std[i] == std_vector[i]);
    }

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

    // TODO make sure we can fix the same template to a different value

    // now release a template
    fit_calc->ReleaseTemplate(fit_calc->GetNTemplates() - 1);
    user_params = Eigen::RowVectorXd::Ones(templates.size() * dims[0]);
    assert((total - fit_calc->U(user_params)).isZero(0));

    user_params(4) = 1.5;

    fit::Minuit2TemplateFitter fitter(fit_calc, 3);
    fitter.SetPrintLevel(0);
    auto result = fitter.Fit(fit_calc->Predict(user_params));

    user_params(4) = 1;

    std::vector<double> p;
    std::vector<double> c;
    std::vector<double> c2;
    std::vector<double> c3;
    std::vector<double> std_params(user_params.size()-1);
    auto eigen_params = user_params;
    double nfits = 5;
    for(auto i = 0u; i < nfits; i++) {
        auto delta = 0. + (20 / nfits) * i;
        fit_calc->FixTemplate(user_params.size()-1, delta);
        auto r = fitter.Fit(total);
        p.push_back(delta);
        c.push_back(r.fun_val);

        for(auto i = 0u; i < std_params.size()-1; i++) {
            std_params[i] = r.params(i);

        }
        std::cout << r.params.transpose() << std::endl;
        std_params[std_params.size()-1] = delta;

        c2.push_back(fit_calc->fun(fit::detail::STDToEigen(std_params), total));
        c3.push_back(fitter(std_params));
    }

    for(auto i = 0u; i < nfits; i++) {
        std::cout << std::setw(20) << p[i] << "|";
    }
    std::cout << std::endl;
    for(auto i = 0u; i < nfits; i++) {
        std::cout << std::setw(20) << c[i] << "|";
    }
    std::cout << std::endl;
    for(auto i = 0u; i < nfits; i++) {
        std::cout << std::setw(20) << c2[i] << "|";
    }
    std::cout << std::endl;
    for(auto i = 0u; i < nfits; i++) {
        std::cout << std::setw(20) << c3[i] << "|";
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------" << std::endl;
    std::cout << result.params.transpose() << std::endl;
    std::cout << result.plus_one_sigma_errors.transpose() << std::endl;
    std::cout << result.minus_one_sigma_errors.transpose() << std::endl;
    std::cout << result.fun_val << std::endl;
    std::cout << result.fun_calls << std::endl;
}


