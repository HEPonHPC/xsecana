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

nominal_loader.Go()
calibration_loader.Go()
lightyeild_up_loader.Go()
lightyeild_down_loader.Go()

"""
analysis contains events DataFrames in memory that have yet to be histogrammed
Have reductions happened? Probably should have
xsecana.OptimizeBinning can run the fit in parallel with MPI
"""
def objective_function(bins):
    histogrammed_analysis = analysis.Histogram(bins)
    return histogrammed_analysis.Result(
        xsecana.SimpleQuadSum
    )[1].sumw2() # optimize against some total error

optimized_bins = xsecana.OptimizeBinning(
    obj_func = objective_function,
)