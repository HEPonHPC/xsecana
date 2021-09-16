import xsecana
import pandana
import numuccinc
import h5py

from functools import partial

loader = pandana.Loader('nominal_file.h5')

def make_unfolder(weight):
    return xsecana.IterativeUnfolder(
        unfolding_matrix=pandana.Spectrum(
            loader,
            VarsToNDVar([numuccinc.kTrueVar,
                         numuccinc.kRecoVar]),
            kSignalCut,
            weight=weight,
        )
    )
def make_unfolder_multiverse(loader, weights):
    return [make_unfolder(loader, weight) for weight in weights]

unfolders_genie_multiverse = xsecana.Systematic(
    [make_unfolder(loader, weight)
     for weight in numuccinc.genie_multiverse_weights]
)
unfolders_flux_multiverse = xsecana.Systematic(
    [make_unfolder(loader, weight)
    for weight in numuccinc.flux_multiverse_weights]
)

output = h5py.File('output.h5')

# reductions happen internally before saved
unfolders_genie_multiverse.SaveTo(output, 'genie')
unfolders_flux_multiverse.SaveTo(output, 'flux')

