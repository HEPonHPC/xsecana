import xsecana
import pandana
from NOvAPandAna.weights import genie_multiverse_weights

from myanalysis import *
from ideal_python_interface_simple_analysis import make_cross_section
"""comments like this are about what the framework is/should/could be doing"""
# other general comments are like this


nominal_loader = pandana.Loader('nominal_file.h5')
calibration_loader = pandana.Loader('calibration_file.h5')
lightyeild_up_loader = pandana.Loader('lightyeild_up_file.h5')
lightyeild_down_loader = pandana.Loader('lightyeild_down_file.h5')

analysis = xsecana.Analysis(
    data = pandana.Spectrum(data_loader, kSelection, kVar),
    nominal = make_cross_section(nominal_loader),
    systematics={
        'calibration_shape': xsecana.Systematic(
            make_cross_section(calibration_loader),
        ),
        'lightyeild': xsecana.Systematic(
            make_cross_section(lightyeild_up_loader),
            make_cross_section(lightyeild_down_loader)
        ),
        'genie_multiverse' : xsecana.Systematic(
            make_cross_section_multiverse(
                nominal_loader,
                genie_multiverse_weights,
            )
        )
    },
)

optimization_var_nominal = pandana.Spectrum(
    nominal_loader,
    kBasicSelection,
    kElectronQuality,
)
optimization_var_calibration_shape = pandana.Spectrum(
    calibration_loader,
    kBasicSelection,
    kElectronQuality,
)
optimization_var_lightyeild = xsecana.Systematic(
    pandana.Spectrum(
        lightyeild_up_loader,
        kBasicSelection,
        kElectronQuality,
    ),
    pandana.Spectrum(
        lightyeild_down_loader,
        kBasicSelection,
        kElectronQuality,
    )
)
optimization_var_genie_multiverse = xsecana.Systematic(
    [pandana.Spectrum(nominal_loader,
                      kLasicSelection,
                      kElectronQuality,
                      weight=weight)
     for weight in genie_multiverse_weights]
)
optimization_var_systematics = {
    'calibration_shape' : optimization_var_calibration_shape,
    'lightyeild' : optimization_var_lightyeild,
    'genie_multiverse' : optimize_var_genie_multiverse,
}
nominal_loader.Go()
calibration_loader.Go()
lightyeild_up_loader.Go()
lightyeild_down_loader.Go()

"""
xsecana.OptimizeCut can run optimization steps in parallel with MPI
"""
def objective_function(nominal_cut : pd.Series,
                       systematics_cut : pd.Series) -> float:
    return xsecana.Analysis(
        nominal=analysis.NominalMeasurement().Cut(nominal_cut),
        systematic=[
            analysis.SystematicMeasurement(name).Cut(systematics_cut[name]) for name in analysis.SystematicNames()
        ]
    ).Histogram(kBinning).Result(xsecana.SimpleQuadSum)[1]

optized_cut = xsecana.OptimizeCut(
    optimization_var_nominal = optimization_var_nominal,
    optimization_var_systematics = optimization_var_systematics,
    range=(0,1),
    nsteps = 100,
    obj_func = objective_function,
)