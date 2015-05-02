from __future__ import print_function, division
import os.path
from functools import partial

from ...external.qt import QtGui
from ...external.qt.QtCore import Qt

from ...qt.widgets.histogram_widget import HistogramWidget, _hash
from ...qt.widgets.mpl_widget import defer_draw
from ...qt.widget_properties import (ButtonProperty, FloatLineProperty,
                                 CurrentComboProperty,
                                 ValueProperty, connect_int_spin)

from ...qt.qtutil import load_ui
from .client import SeqClient

__all__ = ['SeqWidget']


WARN_SLOW = 10000000


class SeqWidget(HistogramWidget):

    LABEL = "Sequence Reader"
    _property_set = HistogramWidget._property_set + 'chrom idn'.split()
    chrom = CurrentComboProperty('ui.chromCombo', 'Chromosome')
    idn = CurrentComboProperty('ui.idnCombo', 'Sample')

    def _load_ui(self):
        self.ui = load_ui(os.path.join(os.path.dirname(__file__), 'SeqWidget.ui'), self.option_widget)

    def _setup_client(self):
        self.client = SeqClient(self._data,
                                self.central_widget.canvas.fig,
                                artist_container=self._container)

    def __str__(self):
        return "Sequence Reader"

    def _connect(self):
        ui = self.ui
        ui.chromCombo.currentIndexChanged.connect(self._update_chrom)
        ui.idnCombo.currentIndexChanged.connect(self._update_idn)
        super(SeqWidget, self)._connect()

    @defer_draw
    def _update_attributes(self):
        """Repopulate the combo box that selects the quantity to plot"""
        combo = self.ui.attributeCombo
        chrom = self.ui.chromCombo
        idn = self.ui.idnCombo
        component = self.component
        new_pos = self.client.component or component
        new_chrom = chrom
        new_idn = idn

        combo.blockSignals(True)
        combo.clear()

        # implementation note:
        # PySide doesn't robustly store python objects with setData
        # use _hash(x) instead
        model = QtGui.QStandardItemModel()
        data_ids = set(_hash(d) for d in self._data)
        self._component_hashes = dict((_hash(c), c) for d in self._data
                                      for c in d.components)

        chr_model = QtGui.QStandardItemModel()
        idn_model = QtGui.QStandardItemModel()

        found_pos = False
        found_chrom = False
        found_idn = False
        for d in self._data:
            if d not in self._container:
                continue
            item = QtGui.QStandardItem(d.label)
            item.setData(_hash(d), role=Qt.UserRole)
            assert item.data(Qt.UserRole) == _hash(d)
            item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            model.appendRow(item)
            for c in d.visible_components:
                if not d.get_component(c).numeric:
                    continue
                if c is new_pos:
                    found_pos = True
                if c.label in ['CHROM', 'chromosome', 'chr', 'CHR',
                               'CHROMOSOME', 'chrom', 'Chromosome', 'Chrom']:
                    for chr in d.get_component(c)._categories:
                        if chr is new_chrom:
                            found_chrom = True
                        itemc = QtGui.QStandardItem(chr)
                        item.setData(_hash(c), role=Qt.UserRole)
                        chr_model.appendRow(itemc)
                if c.label in ['Patient', 'patient', 'PATIENT',
                               'Sample', 'sample', 'SAMPLE']:
                    for ids in d.get_component(c)._categories:
                        if idn is new_idn:
                            found_idn = True
                        itemi = QtGui.QStandardItem(ids)
                        idn_model.appendRow(itemi)
                item = QtGui.QStandardItem(c.label)
                item.setData(_hash(c), role=Qt.UserRole)
                model.appendRow(item)
        combo.setModel(model)
        chrom.setModel(chr_model)
        idn.setModel(idn_model)

        # separators below data items
        for i in range(combo.count()):
            if combo.itemData(i) in data_ids:
                combo.insertSeparator(i + 1)

        combo.blockSignals(False)
        chrom.blockSignals(False)
        idn.blockSignals(False)

        if found_pos:
            self.component = new_pos
        else:
            combo.setCurrentIndex(2)  # skip first data + separator
        self._set_attribute_from_combo()

        if found_chrom:
            self.component = new_chrom
        else:
            combo.setCurrentIndex(1)  # skip first data + separator
        self._update_chrom()

        if found_idn:
            self.component = new_idn
        else:
            combo.setCurrentIndex(1)  # skip first data + separator
        self._update_idn()

    @property
    def chromosome(self):
        combo = self.ui.chromCombo
        index = combo.currentIndex()
        return self._component_hashes.get(combo.itemData(index), None)

    @property
    def sample(self):
        combo = self.ui.idnCombo
        print(combo)
        print(dir(combo))
        index = combo.currentIndex()
        print(index)
        return combo.itemData(index).label

    @defer_draw
    def _set_attribute_from_combo(self, *args):
        self.client.set_component(self.component, 'p')
        self.update_window_title()

    @defer_draw
    def _update_chrom(self, *args):
        sele = self.chrom
        self.client.chrom = sele

    @defer_draw
    def _update_idn(self, *args):
        sele = self.idn
        self.client.idn = sele

    @defer_draw
    def add_data(self, data):
        """ Add data item to combo box.
        If first addition, also update attributes """

        if self.data_present(data):
            return True

        if data.size > WARN_SLOW and not self._confirm_large_data(data):
            return False

        self.client.add_layer(data)
        self._update_attributes()
        self._update_minmax_labels()

        return True