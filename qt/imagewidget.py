from PyQt4.QtCore import Qt, SIGNAL, QVariant
from PyQt4.QtGui import *

import numpy as np

import cloudviz as cv
import cloudviz.message as msg
from cloudviz.image_client import ImageClient

from qtutil import mpl_to_qt4_color, qt4_to_mpl_color
from ui_imagewidget import Ui_ImageWidget
from linker_dialog import LinkerDialog

class ImageWidget(QWidget, cv.HubListener):
    def __init__(self, data):
        QWidget.__init__(self)
        self.ui = Ui_ImageWidget()
        self.ui.setupUi(self)

        self.client = ImageClient(data,
                                  self.ui.mplWidget.canvas.fig,
                                  self.ui.mplWidget.canvas.ax)

        self.connect()
        self.init_widgets()
        self.set_data(0)
        self.client.set_data(self.client.data[0])

    def init_widgets(self):
        self.ui.imageSlider.hide()
        self.ui.sliceComboBox.hide()
        self.ui.sliceComboBox.addItems(["xy", "xz", "yz"])
        for d in self.client.data:
            self.add_data(d)

    def add_data(self, data):
        if len(data.shape) not in [2,3]:
            return
        self.ui.displayDataCombo.addItem(data.label, userData = QVariant(data))

    def set_data(self, index):
        data = self.ui.displayDataCombo.itemData(index).toPyObject()
        self.client.set_data(data)
        self.ui.displayDataCombo.setCurrentIndex(index)
        self.set_attribute_combo(data)
        if not self.client.is_3D:
            self.ui.imageSlider.hide()
            self.ui.sliceComboBox.hide()
            self.ui.orientationLabel.hide()
        else:
            self.ui.imageSlider.show()
            self.ui.sliceComboBox.show()
            self.ui.orientationLabel.show()
        self.set_slider_range()

    def set_attribute(self, index):
        att = self.ui.attributeComboBox.itemText(index)
        att = str(att) # cast from QString
        self.client.set_attribute(att)
        self.ui.attributeComboBox.setCurrentIndex(index)

    def set_attribute_combo(self, data):
        combo = self.ui.attributeComboBox
        combo.currentIndexChanged.disconnect(self.set_attribute)
        combo.clear()
        fields = data.components.keys()
        combo.addItems(fields)
        combo.currentIndexChanged.connect(self.set_attribute)

    def set_slider(self, index):
        self.client.slice_ind = index
        self.ui.imageSlider.setValue(index)

    def set_orientation(self, ori):
        self.client.set_slice_ori(ori)
        self.ui.sliceComboBox.setCurrentIndex(ori)
        self.set_slider_range()

    def set_slider_range(self):
        self.ui.imageSlider.setRange(*self.client.slice_bounds())

    def connect(self):
        ui = self.ui
        client = self.client

        ui.displayDataCombo.currentIndexChanged.connect(self.set_data)
        ui.attributeComboBox.currentIndexChanged.connect(self.set_attribute)
        ui.sliceComboBox.currentIndexChanged.connect(self.set_orientation)

        #connect MPL draw widget mouse events (RMB) to color map manipulation
        ui.imageSlider.sliderMoved.connect(self.set_slider)
        ui.mplWidget.rightDrag.connect(self.set_norm)


    def register_to_hub(self, hub):
        self.client.register_to_hub(hub)
        dc_filt = lambda x:x.sender is self.client._data

        hub.subscribe(self,
                      msg.DataCollectionAddMessage,
                      handler=lambda x:x.add_data(x.sender.data),
                      filter = dc_filt)
        hub.subscribe(self,
                      msg.DataCollectionDeleteMessage,
                      handler=lambda x:x.remove_data(x.sender.data),
                      filter = dc_filt)

    def remove_data(self, data):
        combo = self.ui.displayDataCombo
        for item in range(combo.count()):
            if combo.itemData(item).toPyObject() is data:
                combo.removeItem(item)
                break

    def set_norm(self, x, y):
        if self.client._image is None:
            return
        lo = np.min(self.client._image)
        hi = np.max(self.client._image)
        theta = np.pi - max(min(y, 1), 0) * np.pi
        ra = hi - lo
        bias = lo + ra / 2 * x
        vmin = bias - ra * np.tan(theta)
        vmax = bias + ra * np.tan(theta)
        self.client.set_norm(vmin, vmax)

if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    win = QMainWindow()

    data, subset = cv.example_data.simple_cube()
    d2, s2 = cv.example_data.simple_image()

    dc = cv.DataCollection([data, d2])
    image_client = ImageWidget(dc)

    hub = cv.Hub(data, subset, d2, s2, dc, image_client)

    win.setCentralWidget(image_client)
    win.show()
    #image_client.client.set_norm(1, 2000)
    image_client.client.set_cmap('hot')
    sys.exit(app.exec_())
