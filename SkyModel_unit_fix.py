## This is a very quick fix for a bug that was present before Gammapy v.2.0.1
## Redefine SkyModel init function
from gammapy.utils.scripts import make_name

def __init__2(
    self,
    spectral_model,
    spatial_model=None,
    temporal_model=None,
    name=None,
    apply_irf=None,
    datasets_names=None,
    covariance_data=None,
):
    self.spatial_model = spatial_model
    self.spectral_model = spectral_model
    self.temporal_model = temporal_model
    self._name = make_name(name)

    if apply_irf is None:
        apply_irf = self._apply_irf_default.copy()

    self.apply_irf = apply_irf
    self.datasets_names = datasets_names
    #self._check_unit()

    gammapy.modeling.models.core.ModelBase.__init__(self=self,covariance_data=covariance_data)
SkyModel.__init__=__init__2