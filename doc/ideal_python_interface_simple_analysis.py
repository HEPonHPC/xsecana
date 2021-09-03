import xsecana
import pandana

# user-defined function
def make_cross_section(loader):
    efficiency = xsecana.SimpleEfficiency(
        numerator = pandana.Spectrum(loader, kVar, kSelection & kReconstructedSignal),
        denominator = pandana.Spectrum(loader, kVar, kAllSignal)
    )

    flux = xsecana.SimpleFlux(
        pandana.Spectrum(loader, kNeutrinoEnergy, kNuebarCC)
    )

    unfolder = xsecana.IterativeUnfolder(
        matrix=pandana.Spectrum(
            loader,
            kRecoVsTrue,
            kReconstructedSignal
        )
    )

    signal_estimator = xsecana.SimpleSignalEstimator(
        signal = pandana.Spectrum(loader, kVar, kSelection & kReconstructedSignal),
        background = pandana.Spectrum(loader, kVar, kSelection & kReconstructedBackground)
    )

    cross_section = xsecana.CrossSection (
        efficiency = efficiency,
        flux = flux,
        unfolder = unfolder,
        signal_estimator = signal_estimator,
        ntargets = 1e4,
    )

    return cross_section

nominal_loader = pandana.Loader('nominal_file.h5')
calibration_loader = pandana.Loader('calibration_file.h5')
lightyeild_loader = pandana.Loader('lightyeild_file.h5')

xsecana.Analysis(
    nominal = make_cross_section(nominal_loader),
    systematics = {
        'calibration': make_cross_section(calibration_loader),
        'lightyeild': make_cross_section(lightyeild_loader)
    }
    uncertainty_propagator = xsecana.SimpleQuadSum,
)

nominal_loader.Go()
calibration_loader.Go()
lightyeild_loader.Go()

# should be able to save data to file
# so analysis could be performed as a separate step from the selection.
# Useful for current analysis workflow where a grid job creates many files
# that are aggregated after the selection
xsecana.SaveTo('analysis_file.h5')

# but also we should be able to just use it right away (HPC + MPI jobs)
xsecana.MPIReduce(MPI.COMM_WORLD)
central_value_result = xsecana.Result()
result_uncertainty = xsecana.TotalFractionalUncertainty()

# save results for plotting later
