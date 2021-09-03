import xsecana
import pandana

import myanalysis

loader = pandana.Loader('file.h5')

# user-defined function that will interact with framework
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

# maybe we want to save to file first and do the study later
signal_estimator.SaveTo('signal_estimator.h5')

# or maybe we want to do it right away

# is it possible to somehow do MPI reductions with some generality?
# signal_estimator.Reduce(MPI.COMM_WORLD)

signal_estimator.DoMyInterestingStudy()


