import xsecana
import pandana

# user library
from myanalysis import *

"""comments like this are about what the framework is/should be doing"""
# other general comments are like this

# user-defined function
def make_cross_section(loader):
    # pre-packaged from framework
    """interface with pandana"""
    efficiency = xsecana.SimpleEfficiency(
        numerator = pandana.Spectrum(loader, kVar, kSelection & kReconstructedSignal),
        denominator = pandana.Spectrum(loader, kVar, kAllSignal)
    )

    # pre-packaged from framework
    """interface with pandana"""
    flux = xsecana.SimpleFlux(
        pandana.Spectrum(loader, kNeutrinoEnergy, kNuebarCC)
    )

    # pre-packaged from framework
    """interface with pandana"""
    unfolder = xsecana.IterativeUnfolder(
        matrix=pandana.Spectrum(
            loader,
            kRecoVsTrue,
            kReconstructedSignal
        )
    )

    # user-defined function that will interact with framework
    """
    interface with pandana
    User's SignalEstimator inherits from xsecana.ISignalEstimator
    """
    signal_estimator = myanalysis.SignalEstimator(
        templates={
            'signal': pandana.Spectrum(
                loader, kVar, kSelection & kReconstructedSignal
            ),
            'background_1': pandana.Spectrum(
                loader, kVar, kSelection & kReconstructedBackground_1
            ),
            'background_2': pandana.Spectrum(
                loader, kVar, kSelection & kReconstructedBackground_2
            ),
        }

    )

    """
    - Given these components, this calculates the cross section and
      keeps everything together.
    - A TemplateFitEstimator actually depends on all systematic samples
      to build a covariance matrix.
    - 1-sample of a CrossSection with a TemplateFitSignalEstimator is incomplete,
      it cannot calculate a cross section
      - Keeping all of these as a CrossSection object is associating all of the 
        right components though.
      - Is still useful for correlating the flux and efficiency uncertainties for 
        each sample
      - CrossSection objects should not assume they are complete
    """
    cross_section = xsecana.CrossSection(
        efficiency = efficiency,
        flux = flux,
        unfolder = unfolder,
        signal_estimator = signal_estimator,
        ntargets = 1e4,
    )

    return cross_section

data_loader = pandana.Loader('data_file.h5')
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
"""
data = pandana.Spectrum(data_loader, kSelection, kVar)
nominal = make_cross_section(nominal_loader)
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

nominal_loader.Go()
calibration_loader.Go()
lightyeild_up_loader.Go()
lightyeild_down_loader.Go()

# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
"""
- MPI reductions happen here
- Can make framework save all of the user's objects, but
  how do we load them later? If python, we can just pickle them
- Events can still be saved in pd.DataFrames

xsecana.SaveAnalysis(
  group: h5py.Group,
  nominal : xsecana.IMeasurement,
  systematics : dict(str, xsecana.Systematic),
  data : xsecana.Hist,
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

myanalysis.LoadMyMeasurement(
  handle: h5py.File|h5py.Group|h5py.Dataset,
  name : str
)-> xsecana.IMeasurement

xsec
"""
# if starting from analysis file on disk
with h5py.File('analysis_file.h5', 'r') as f:
    data, nominal, systematics = xsecana.LoadAnalysis(
        f.get('myanalysis'),
        myanalysis.LoadMyMeasurement,
    )


# but also we should be able to just use it right away (HPC + MPI jobs)

# selected events are aggregated behind the scenes when they're needed

"""
xsecana.TemplateFit.TotalUncertainty(
  binning : np.array,
  nominal : xsecana.IMeausurement,
  systematics : dict(str, xsecana.Systematic),
  data : xsecana.Array,
)-> (xsecana.Hist, xsecana.Systematic)
"""
central_value, uncertainty = xsecana.TemplateFit.TotalUncertainty(
    myanalysis.TemplateBinning,
    nominal,
    systematics,
    data
)
