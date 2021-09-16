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
    }
)

nominal_loader.Go()
calibration_loader.Go()
lightyeild_up_loader.Go()
lightyeild_down_loader.Go()

# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
"""
- MPI reductions happen here
- Can make framework save all of the user's objects, but
  how do we load them later?
"""
analysis.SaveTo(h5py.File('analysis_file.h5'),
                'myanalysis')

# if starting from analysis file on disk
analysis = xsecana.Analysis.LoadFrom(
    myanalysis.LoadMyMeasurement,
    h5py.File('analysis_file.h5', 'r'),
    'myanalysis',
)
"""
A callable like myanlaysis.LoadMyMeasurement is how analyzers tell the framework
how to load their classes. It is required to have a function signature like:
(handle: hdf5.File|hdf5.Group|hdf5.Dataset, group_name: str)-> xsecana.IMeasurement 
"""


# but also we should be able to just use it right away (HPC + MPI jobs)

# selected events are aggregated behind the scenes when they're needed
"""
If we hadn't just saved the results, MPI reductions would happen here instead
Analysis object calls the TemplateFitUncertainty function with nominal, systematics, and data
This is where the events can be optionally histogrammed

expecting a callable like:
fun(
  nominal     : xsecana.IMeasurement, 
  systematics : dict(str, xsecana.IMeasurement),
  data        : xsecana.Hist
) -> (xsecana.Hist, xsecana.Systematic)  
"""
central_value, uncertainty = analysis.Result(myanalysis.TemplateFitUncertainty) # or xsecana.TemplateFitUncertainty

"""
What would it look like if there just wasn't an analysis object? It doesn't do a a whole lot
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
Save to file
"""
nominal.SaveTo(output, 'nominal')
data.SaveTo(output, 'data')
systematics.SaveTo('systematics')

central_value, uncertainty = myanalysis.TemplateFitUncertainty(
    nominal,
    systematics,
    data
)

"""
So Analysis is saving two lines of code in this case. 
"""

