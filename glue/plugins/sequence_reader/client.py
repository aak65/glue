from __future__ import print_function, division
from abc import ABCMeta, abstractproperty, abstractmethod
import numpy as np
from functools import partial

from ...external import six
from ...core.data import IncompatibleAttribute, CategoricalComponent
from ...core.subset import RangeSubsetState
from ...core.edit_subset_mode import EditSubsetMode
from ...core.callback_property import (CallbackProperty, add_callback,
                                      delay_callback)

from ...clients.histogram_client import HistogramClient, UpdateProperty
from ...clients.layer_artist import (ChangedTrigger, HistogramLayerArtist)


__all__ = ['SeqClient', 'SeqLayerBase', 'SeqLayerArtist']


@six.add_metaclass(ABCMeta)
class SeqLayerBase(object):
    lo = abstractproperty()     # lo-cutoff for bin counting
    hi = abstractproperty()     # hi-cutoff for bin counting
    nbins = abstractproperty()  # number of bins
    xlog = abstractproperty()   # whether to space bins logarithmically
    chrom = abstractproperty()  # chromosome to display
    idn = abstractproperty()     # for use with the categorical ROI

    @abstractmethod
    def get_data(self):
        """
        Return array of bin counts
        """
        pass


class SeqLayerArtist(HistogramLayerArtist, SeqLayerBase):
    _property_set = HistogramLayerArtist._property_set + 'chrom idn'.split()
    chrom = ChangedTrigger()
    idn = ChangedTrigger()

    def _calculate_histogram(self):
        """Recalculate the histogram, creating new patches"""
        self.clear()
        try:
            data = self.layer[self.att].ravel()
            chr = self.layer[self.chrom].ravel()
            idn = self.layer[self.idn].ravel()
            if not np.isfinite(data).any():
                return False
        except IncompatibleAttribute as exc:
            self.disable_invalid_attributes(*exc.args)
            return False

        if data.size == 0:
            return

        if self.lo > np.nanmax(data) or self.hi < np.nanmin(data):
            return
        if self.xlog:
            data = np.log10(data)
            rng = [np.log10(self.lo), np.log10(self.hi)]
        else:
            rng = self.lo, self.hi

        # to_patch = [x for i, x in enumerate(data) if chr[i] == 1 and idn[i] == 1]
        nbinpatch = self._axes.hist(to_patch,
                                    bins=self.nbins,
                                    range=rng)
        self._y, self.x, self.artists = nbinpatch
        return True


class SeqClient(HistogramClient):

    """
    A client class that uses matplotlib to visualize tables as scatter plots.
    """
    chrom = CallbackProperty()
    idn = CallbackProperty()
    layer_artist_class = SeqLayerArtist

    def _connect(self):
        add_callback(self, 'chrom', partial(self.set_component, self._component, 'c'))
        add_callback(self, 'idn', partial(self.set_component, self._component, 'i'))

    def set_component(self, component, coord, *args):
        """
        Redefine which component gets plotted

        Parameters
        ----------
        component: ComponentID
            The new component to plot
        """
        if self._component is component:
            return

        iscat = lambda x: isinstance(x, CategoricalComponent)

        def comp_obj():
            # the current Component (not ComponentID) object
            x = list(self._get_data_components('x'))
            if x:
                x = x[0]
            return x

        prev = comp_obj()
        old = self.nbins

        first_add = self._component is None
        self._component = component
        cur = comp_obj()

        if first_add or iscat(cur):
            self._auto_nbin()

        # save old bins if switch from non-category to category
        if prev and not iscat(prev) and iscat(cur):
            self._saved_nbins = old

        # restore old bins if switch from category to non-category
        if iscat(prev) and cur and not iscat(cur) and self._saved_nbins is not None:
            self.nbins = self._saved_nbins
            self._saved_nbins = None

        self.sync_all()
        self._relim()

    def apply_roi(self, roi):
        x, _ = roi.to_polygon()
        lo = min(x)
        hi = max(x)

        # expand roi to match bin edges
        bins = self.bins

        if lo >= bins.min():
            lo = bins[bins <= lo].max()
        if hi <= bins.max():
            hi = bins[bins >= hi].min()

        if self.xlog:
            lo = 10 ** lo
            hi = 10 ** hi

        state = RangeSubsetState(lo, hi)
        state.att = self.component
        mode = EditSubsetMode()
        visible = [d for d in self.data if self.is_layer_visible(d)]
        focus = visible[0] if len(visible) > 0 else None
        mode.update(self.data, state, focus_data=focus)



