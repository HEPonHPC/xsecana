import xsecana
import pandana
from NOvAPandAna.weights import genie_multiverse_weights

from myanalysis import *

"""comments like this are about what the framework is/should/could be doing"""
# other general comments are like this

# user-defined function
def make_cross_sections(loader, weights: List(pandana.Weight)=None):
    if type(weights) is not list: weights = [weights]

    selected_signal = [
        pandana.Spectrum(
            loader,
            kVar,
            kSelection & kReconstructedSignal,
            weight=weight
        )
        for weight in weights
    ]
    selected_background = [
        pandana.Spectrum(
            loader,
            kVar,
            kSelection & kReconstructedBackground,
            weight=weight
        )
        for weight in weights
    ]
    all_signal = [
        pandana.Spectrum(
            loader,
            kVar,
            kAllSignal,
            weight=weight
        )
        for weight in weights
    ]

    unfolding_matrix = [
        pandana.Spectrum(
            loader,
            kRecoVsTrue,
            kReconstructedSignal,
            weight=weight,
        )
        for weight in weights
    ]

    flux = NOvAPandAna.DeriveFlux(pdg=-14)
    # fill spectra
    loader.Go()

    """
    - Given these components, this calculates the cross section and
      keeps everything together.
    - Together, these components can calculate a cross section for this 
      particular sample independently. It is complete.
    """
    cross_sections = [xsecana.CrossSectionCalculator(
        efficiency=xsecana.SimpleEfficiencyCalculator(
            selected_signal[iweight]._df,
            all_signal[weight]._df
        ),
        flux=xsecana.SimpleFluxCalculator(
            flux._df
        ),  # xsecana.IFlux
        unfolder=xsecana.IterativeUnfolder(
            unfolding_matrix[weight]._df
        ),
        signal_estimator=xsecana.SimpleSignalEstimationCalculator(
            selected_background[iweight]._df
        ),
        ntargets=kNTargets,  # int
    )
        for iweight in range(len(weights))
    ]
    return cross_sections[0] if len(weights) == 1 else cross_sections

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

data_spectrum = pandana.Spectrum(data_loader, kSelection, kVar)
data_loader.Go()

data = xsecana.Array(data_spectrum._df)
nominal = make_cross_section(nominal_loader)
systematics = {
    'calibration_shape': xsecana.Systematic(
        make_cross_sections(calibration_loader),
    ),
    'lightyeild': xsecana.Systematic(
        make_cross_sections(lightyeild_up_loader),
        make_cross_sections(lightyeild_down_loader)
    ),
    'genie_multiverse': xsecana.Systematic(
        make_cross_sections(
            nominal_loader,
            weights = genie_multiverse_weights,
        )
    )
}

# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
"""
- MPI reductions happen on save
xsecana.SaveAnalysis(
  group : hdf5.Group,
  data : llama.array,
  nominal : xsecana.IMeasurement,
  systematics : dict(str, xsecana.Systematic),
)->None
"""
with h5py.File('analysis_file.h5', 'w') as f:
    xsecana.SaveAnalysis(
        f.create_group('myanalysis'),
        data,
        nominal,
        systematics,
    )

"""
A callable like myanlaysis.LoadMyMeasurement is how analyzers tell the framework
how to load their classes. It is required to have a function signature like:
(handle: hdf5.File|hdf5.Group|hdf5.Dataset, group_name: str)-> xsecana.IMeasurement 
"""
# if starting from analysis file on disk
with h5py.File('analysis_file.h5', 'r') as f:
    data, nominal, systematics = xsecana.LoadAnalysis(
        f.get('myanalysis'),
        myanalysis.LoadMyMeasurement,
    )

# but also we should be able to just use it right away (HPC + MPI jobs)

"""
Events are aggregated behind the scenes when they're needed
for these calculations
"""
# maybe we want an individual uncertainty when events are summarized in
# some given binning
central_value, uncertainty = xsecana.SimpleQuadSum.AbsoluteUncertainty(
    data,
    nominal,
    systematics['calibration_shape'],
    myanalysis.kBinning,
)

# or maybe we want the total uncertainty on final result when events are
# summarized in some given binning
central_value, uncertainty = xsecana.SimpleQuadSum.TotalAbsoluteUncertainty(
    data,
    nominal,
    systematics,
    myanalysis.kBinning,
)

