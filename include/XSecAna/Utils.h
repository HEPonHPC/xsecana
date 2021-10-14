#pragma once

#include <Eigen/Dense>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "XSecAna/Type.h"

namespace xsec {
    typedef Eigen::ArrayXd Array;
    typedef Eigen::ArrayXXd Array2D;
    typedef Eigen::Tensor<double, 3> Array3D;

    typedef Eigen::Ref<Array> ArrayRef;
    typedef Eigen::Map<const Array> ArrayMap;
    namespace root {
        inline std::string MakeUnique(const std::string & base) {
            static int N = 0;
            return base + std::to_string(N++);
        }

        inline std::unique_ptr<TH1> LoadTH1(TDirectory * dir, const std::string & name) {
            auto ret = std::unique_ptr<TH1>((TH1 *) dir->Get(name.c_str()));
            ret->SetDirectory(0); // Tell dir it doesn't own this object anymore
            // or else it gets deleted when the file is closed...
            return ret;
        };

        struct TH1Props {
            TH1Props() = default;

            TH1Props(const TH1 * h, const char * _name = "")
                    : axes({h->GetXaxis(),
                            h->GetYaxis(),
                            h->GetZaxis()}),
                      dims(h->GetDimension()),
                      name(_name),
                      entries(h->GetEntries()) {
                if (dims == 1) {
                    nbins_and_uof = h->GetNbinsX() + 2;
                } else if (dims == 2) {
                    nbins_and_uof = (h->GetNbinsX() + 2) *
                                    (h->GetNbinsY() + 2);
                } else {
                    nbins_and_uof = (h->GetNbinsX() + 2) *
                                    (h->GetNbinsY() + 2) *
                                    (h->GetNbinsZ() + 2);
                }
            }


            std::vector<const TAxis *> axes;
            int dims;
            const char * name;
            unsigned int nbins_and_uof;
            unsigned int entries;
        };

        inline ArrayMap MapContentsToEigen(const TH1 * h) {
            if (h->GetDimension() == 1) {
                return ArrayMap(dynamic_cast<const TH1D *>(h)->GetArray(),
                                h->GetNbinsX() + 2);
            } else if (h->GetDimension() == 2) {
                return ArrayMap(dynamic_cast<const TH2D *>(h)->GetArray(),
                                (h->GetNbinsX() + 2) *
                                (h->GetNbinsY() + 2));
            } else {
                return ArrayMap(dynamic_cast<const TH3D *>(h)->GetArray(),
                                (h->GetNbinsX() + 2) *
                                (h->GetNbinsY() + 2) *
                                (h->GetNbinsZ() + 2));
            }
        }

        inline Array MapErrorsToEigen(const TH1 * h) {
            if (h->GetDimension() == 1) {
                Array errors(h->GetNbinsX() + 2);
                for (auto i = 0u; i < errors.size(); i++) {
 		     errors(i) = h->GetBinError(i);
                }
                return ArrayMap(errors.data(),
                                errors.size());


            } else if (h->GetDimension() == 2) {
                Array2D errors((h->GetNbinsX() + 2),
                               (h->GetNbinsY() + 2));
                for (auto i = 0u; i < errors.rows(); i++) {
                    for (auto j = 0u; j < errors.cols(); j++) {
                        errors(i, j) = h->GetBinError(i, j);
                    }
                }
                return ArrayMap(errors.data(), errors.size());


            } else {
                Array3D errors(h->GetNbinsX() + 2,
                               h->GetNbinsY() + 2,
                               h->GetNbinsZ() + 2);
                for (auto i = 0u; i <= h->GetNbinsX() + 1; i++) {
                    for (auto j = 0u; j <= h->GetNbinsY() + 1; j++) {
                        for (auto k = 0u; k < h->GetNbinsZ() + 1; k++) {
                            errors(i, j, k) = h->GetBinError(i, j, k);
                        }
                    }
                }
                return ArrayMap(errors.data(),
                                errors.size());
            }
        }

        inline TH1 * ToROOT(const Array & data,
                            TH1Props props) {
            TH1 * h;
            if (props.dims == 1) {
                h = new TH1D(props.name,
                             "",
                             props.axes[0]->GetNbins(),
                             props.axes[0]->GetXbins()->GetArray());
            } else if (props.dims == 2) {
                h = new TH2D(props.name,
                             "",
                             props.axes[0]->GetNbins(),
                             props.axes[0]->GetXbins()->GetArray(),
                             props.axes[1]->GetNbins(),
                             props.axes[1]->GetXbins()->GetArray());
            } else if (props.dims == 3) {
                h = new TH3D(props.name,
                             "",
                             props.axes[0]->GetNbins(),
                             props.axes[0]->GetXbins()->GetArray(),
                             props.axes[1]->GetNbins(),
                             props.axes[1]->GetXbins()->GetArray(),
                             props.axes[2]->GetNbins(),
                             props.axes[2]->GetXbins()->GetArray());
            } else {
                return 0;
            }
            h->SetContent(data.data());
            h->SetEntries(props.entries);
            h->Sumw2();
            return h;
        }

        inline TH1 * ToROOT(const Array & data,
                            const Array & error,
                            TH1Props props) {
            auto ret = ToROOT(data, props);
            ret->SetError(error.data());
            return ret;
        }

        inline TH1 * ToROOTLike(const TH1 * h,
                                const Array & data,
                                const Array & error) {
            auto ret = ToROOT(data, TH1Props(h, ""));
            ret->SetError(error.data());
            return ret;
        }

        inline TH1 * ToROOTLike(const TH1 * h,
                                const Array & data) {
            auto ret = ToROOT(data, TH1Props(h, ""));
            return ret;
        }


        inline void MapToEigen(const TH1 * h, ArrayRef contents, ArrayRef errors) {
            contents = MapContentsToEigen(h);
            errors = MapErrorsToEigen(h);
        }

        inline Array MapBinWidthsToEigen(const TH1Props & props) {
            if (props.dims == 1) {
                Array _bw = Array::Ones(props.axes[0]->GetNbins()+2);
                for (auto i = 1u; i <= props.axes[0]->GetNbins(); i++) {
                    _bw(i) = props.axes[0]->GetBinWidth(i);
                }

                return _bw;
            } else if (props.dims == 2) {
                Array2D _bw = Array2D::Ones(props.axes[0]->GetNbins()+2,
                                            props.axes[1]->GetNbins()+2);
                for (auto i = 1u; i <= props.axes[0]->GetNbins(); i++) {
                    for (auto j = 1u; j <= props.axes[1]->GetNbins(); j++) {
                        _bw(i, j) = props.axes[0]->GetBinWidth(i) *
                                    props.axes[1]->GetBinWidth(j);
                    }
                }
                return ArrayMap(_bw.data(), _bw.size());
            } else {
                Array3D _bw = Eigen::TensorMap<const Array3D>(Array::Ones(props.nbins_and_uof).eval().data(),
                                                              props.axes[0]->GetNbins() + 2,
                                                              props.axes[1]->GetNbins() + 2,
                                                              props.axes[2]->GetNbins() + 2);

                for (auto i = 1u; i <= props.axes[0]->GetNbins(); i++) {
                    for (auto j = 1u; j <= props.axes[1]->GetNbins(); j++) {
                        for (auto k = 1u; k <= props.axes[2]->GetNbins(); k++) {
                            _bw(i - 1, j - 1, k - 1) = props.axes[0]->GetBinWidth(i) *
                                                       props.axes[1]->GetBinWidth(j) *
                                                       props.axes[2]->GetBinWidth(k);
                        }
                    }
                }
                return ArrayMap(_bw.data(), _bw.size());
            }
        }

        inline Array MapBinWidthsToEigen(const TH1 * h) {
            return MapBinWidthsToEigen(TH1Props(h, ""));
        }

    }
}
