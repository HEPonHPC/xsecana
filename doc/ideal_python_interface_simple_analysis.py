import xsecana
import pandana
from NOvAPandAna.weights import genie_multiverse_weights

from myanalysis import *

"""comments like this are about what the framework is/should/could be doing"""
# other general comments are like this

# user-defined function
def make_cross_section(loader, weight=None):
    """interface with pandana"""
    efficiency = xsecana.SimpleEfficiency(
        numerator = pandana.Spectrum(loader, kVar, kSelection & kReconstructedSignal, weight=weight),
        denominator = pandana.Spectrum(loader, kVar, kAllSignal, weight=weight),
        binning=kBinning, # or xsecana.kUnbinned/None?
    )

    """interface with pandana"""
    flux = xsecana.SimpleFlux(
        pandana.Spectrum(loader, kNeutrinoEnergy, kNuebarCC, weight=weight),
        integrated = True,
        binning = kBinning, # or xsecana.kUnbinned/None?
    )

    """interface with pandana"""
    unfolder = xsecana.IterativeUnfolder(
        matrix=pandana.Spectrum(
            loader,
            kRecoVsTrue,
            kReconstructedSignal,
            weight=weight,
        ),
        binning = kBinning, # or xsecana.kUnbinned/None?
    )

    """interface with pandana"""
    signal_estimator = xsecana.SimpleSignalEstimator(
        signal = pandana.Spectrum(loader, kVar, kSelection & kReconstructedSignal, weight=weight),
        background = pandana.Spectrum(loader, kVar, kSelection & kReconstructedBackground, weight=weight),
        binning = kBinning, # or xsecana.kUnbinned/None?
    )
    """
    - Given these components, this calculates the cross section and
      keeps everything together.
    - Together, these components can calculate a cross section for this 
      particular sample independently. It is complete.
    """
    cross_section = xsecana.CrossSection(
        efficiency = efficiency,
        flux = flux, # xsecana.IFlux
        unfolder = unfolder,
        signal_estimator = signal_estimator,
        ntargets = kNTargets, # int
    )
    return cross_section

nominal_loader = pandana.Loader('nominal_file.h5')
calibration_loader = pandana.Loader('calibration_file.h5')
lightyeild_up_loader = pandana.Loader('lightyeild_up_file.h5')
lightyeild_down_loader = pandana.Loader('lightyeild_down_file.h5')

"""
Systematic objects let the user handle multiple associated samples as one
object with the same interface for the three types of systematics
  1. 1-sided
  2. 2-sided
  3. Multiverse
Analysis object keeps a container of Systematics, nominal measurement, and "data"
Can it also propagate binning throughout the entire analysis?
"""
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
When using pandana, the events are held in pd.DataFrames
When and how does optional histogramming happen?
"""


# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
"""
- MPI reductions happen here
- Can make framework save all of the user's objects, but
  how do we load them later? If python, we can just pickle them
- Events can still be saved in pd.DataFrames
"""
analysis.SaveTo(h5py.File('analysis_file.h5'),
                'myanalysis')

# if starting from analysis file on disk
analysis = xsecana.Analysis.LoadFrom(
    myanalysis.LoadMyMeasurement,
    h5py.File('analysis_file.h5', 'r'),
    'myanalysis',
)

# but also we should be able to just use it right away (HPC + MPI jobs)

# selected events are aggregated behind the scenes when they're needed
"""
If we hadn't just saved the results, MPI reductions would happen here instead
"""
histogrammed_analysis = analysis.Histogram(kBinning)
"""
Notes on histogramming:
- Efficiency, SignalEstimator, and Unfolder will be binned similarly
- Flux is unique.
  - If measurement is flux-integrated, the binning will 
    be the same, but the contents of that histogram will all be equal to the 
    flux integral.
  - Otherwise, the flux could be measured in the same space as the analysis
"""

"""
Analysis object calls the SimpleQuadSum uncertainty function with nominal, systematics, and data
expecting a callable like:
fun(
  nominal     : xsecana.IMeasurement, 
  systematics : dict(str, xsecana.Systematic),
  data        : xsecana.Hist
) -> (xsecana.Hist, xsecana.Systematic)
"""
central_value, uncertainty = analysis.Result(xsecana.SimpleQuadSum.TotalAbsoluteUncertainty)



"""
What would it look like if there just wasn't an analysis object?
It doesn't do a a whole lot here.
Remember Systematics are dealing with CrossSection objects here
nominal is also a CrossSection object 
"""
data = pandana.Spectrum(data_loader, kSelection, kVar),
nominal = make_cross_section(nominal_loader),
systematics = {
    'calibration_shape': xsecana.Systematic(
        make_cross_section(calibration_loader),
    ),
    'lightyeild': xsecana.Systematic(
        make_cross_section(lightyeild_up_loader),
        make_cross_section(lightyeild_down_loader)
    ),
    'genie_multiverse': xsecana.Systematic(
        make_cross_section_multiverse(
            nominal_loader,
            genie_multiverse_weights,
        )
    )
}

"""
Save to file. CrossSection object can save user's components, but
doesn't know how to load them, but we would like to be able to load. Same issue here.
"""
nominal.SaveTo(output, 'nominal')
data.SaveTo(output, 'data')
systematics.SaveTo('systematics')

central_value, uncertainty = xsecana.SimpleQuadSum.TotalAbsoluteUncertainty(
    nominal,
    systematics,
    data
)

"""
So Analysis is saving two lines of code in this case. 
"""