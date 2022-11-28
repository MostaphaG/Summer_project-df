
import df3_for_gui as df3
import numpy as np


import os
os.environ['ETS_TOOLKIT'] = 'qt4'
# By default, the PySide binding will be used. If you want the PyQt bindings
# to be used, you need to set the QT_API environment variable to 'pyqt'
#os.environ['QT_API'] = 'pyqt'

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
from pyface.qt import QtGui, QtCore
# Alternatively, you can bypass this line, but you need to make sure that
# the following lines are executed before the import of PyQT:
#   import sip
#   sip.setapi('QString', 2)

from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor


################################################################################
#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.

        # We can do normal mlab calls on the embedded scene.

        grid = np.linspace(-1,1,3)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        fx = xg/np.sqrt(xg**2+yg**2 - zg**2)
        fy = yg/np.sqrt(xg**2+yg**2- zg**2)
        fz = yg/np.sqrt(xg**2+yg**2- zg**2)

        form_1 = df3.form_1_3d(xg, yg, zg, fx, fy, fz)
        stck_coords, red_balls_data, axes_limits = form_1.plot()

        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)

        self.scene.add_actor(stck_coords[0])
        self.scene.add_actor(stck_coords[1])
        self.scene.add_actor(stck_coords[2])
        self.scene.add_actor(stck_coords[3])
        self.scene.add_actor(stck_coords[4])
        self.scene.add_actor(stck_coords[5])

        self.scene.mlab.points3d(red_balls_data[0],
                                 red_balls_data[1],
                                 red_balls_data[2], color = (1,0,0),scale_factor=red_balls_data[3], resolution=36)

        self.scene.mlab.axes(color = (0,0,0), nb_labels = 5, extent = axes_limits, line_width=3.0)
        



    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True # We need this to resize with the parent widget
                )


################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        #QtCore.pyqtRemoveInputHook()
        #import pdb ; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)


if __name__ == "__main__":
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    container.setWindowTitle("Embedding Mayavi in a PyQt4 Application")
    # define a "complex" layout to test the behaviour
    layout = QtGui.QGridLayout(container)

    # put some stuff around mayavi
    label_list = []
    for i in range(3):
        for j in range(3):
            if (i==1) and (j==1):continue
            label = QtGui.QLabel(container)
            label.setText("Your QWidget at (%d, %d)" % (i,j))
            label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
            layout.addWidget(label, i, j)
            label_list.append(label)
    mayavi_widget = MayaviQWidget(container)

    layout.addWidget(mayavi_widget, 1, 1)
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()

    # Start the main event loop.
    app.exec_()