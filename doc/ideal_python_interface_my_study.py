import xsecana
import pandana

import myanalysis

def make_signal_estimator(loader):
    return myanalysis.PurityEstimator(
        signal=pandana.Spectrum(
            loader, kVar, kSelection & kReconstructedSignal
        ),
        background=pandana.Spectrum(
            loader, kVar, kSelection & kReconstructedBackground
        ),
    )

nominal_loader = pandana.Loader('nominal_file.h5')
calibration_loader = pandana.Loader('calibration_file.h5')
lightyeild_loader = pandana.Loader('lightyeild_file.h5')

# user-defined function that will interact with framework
nominal_signal_estimator = make_signal_estimator(nominal_loader)
calibration_signal_estimator = make_signal_estimator(calibration_loader)
lightyeild_signal_estimator = make_signal_estimator(lightyeild_loader)

nominal_loader.Go()
calibration_loader.Go()
lightyeild_loader.Go()

# maybe we want to save to file first and do the study later
signal_estimator.SaveTo('signal_estimator.h5')

# or maybe we want to do it right away

# is it possible to somehow do mpi reductions with some generality?
# signal_estimator.reduce(mpi.comm_world)

lightyeild_uncertainty = xsecana.SimpleQuadSum(
    nominal_signal_estimator,
    lightyeild_signal_estimator,
)

calibration_uncertainty = xsecana.SimpleQuadSum(
    nominal_signal_estimator,
    calibration_signal_estimator,
)



