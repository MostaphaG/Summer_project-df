import os
os.environ['ETS_TOOLKIT'] = 'qt4'

 
from mayavi import mlab
from tvtk.api import tvtk
from tvtk.common import configure_input_data
from tvtk.pyface.api import Scene
from numpy import arange, nonzero, float32, min, max, median, copy, random, shape
from numpy.core.numeric import ravel
 
from traits.api import HasTraits, Instance, on_trait_change, \
    Int, Dict
from traitsui.api import View, Item, VGroup, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor


from pyface.qt import QtGui, QtCore

import sys

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import numpy as np
import df3_for_gui as df3





lft = -5
rght = 5
pts = 10



class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())
    view = View(Item('scene', 
                     editor=SceneEditor(scene_class=MayaviScene),
                     height=250, 
                     width=300, 
                     show_label=False),
                resizable=True
                )
 
    needUpdate = None


    def clear(self):
        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()


 
    def takePlotParametresF1(self, stck_coords, red_balls_data_f1, axes_limits_f1):

        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()
        self.stck_coords = stck_coords
        self.red_balls_data_f1 = red_balls_data_f1
        self.axes_limits_f1 = axes_limits_f1

        
        self.needUpdate = True
 
        self.update_plot_f1()

    def takePlotParametresVF(self, fld_comps, red_balls_data_vf, axes_limits_vf):

        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()
        self.fld_comps = fld_comps
        self.red_balls_data_vf = red_balls_data_vf
        self.axes_limits_vf = axes_limits_vf

        
        self.needUpdate = True
 
        self.update_plot_vf()

    def takePlotParametresF0(self, sc_fld, axes_limits_f0):

        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()
        self.sc_fld = sc_fld
        self.axes_limits_f0 = axes_limits_f0

        
        self.needUpdate = True
 
        self.update_plot_f0()

    def takePlotParametresF2(self, xg, yg, zg, fx, fy, fz, fx_eqn, fy_eqn, fz_eqn):

        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()
        self.x = xg
        self.y = yg
        self.z = zg
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.fx_eqn = fx_eqn
        self.fy_eqn = fy_eqn
        self.fz_eqn = fz_eqn

        
        self.needUpdate = True
 
        self.update_plot_f2()

    def takePlotParametresF3(self, xg3, yg3, zg3, f3, f3_eqn):

        self.scene.mlab.clf()
        self.scene.renderer.remove_all_view_props()
        self.x3 = xg3
        self.y3 = yg3
        self.z3 = zg3
        self.f3 = f3
        self.f3_eqn = f3_eqn


        
        self.needUpdate = True
 
        self.update_plot_f3()

        
 
    @on_trait_change('scene.activated')
    def update_plot_f1(self):

        if not self.needUpdate:
            vtext = tvtk.VectorText()
            vtext.text = 'DFormPy 3D'
            text_mapper = tvtk.PolyDataMapper()
            configure_input_data(text_mapper, vtext.get_output())
            vtext.update()
            p2 = tvtk.Property(color=(0, 0.3, 0.3))
            text_actor = tvtk.Follower(mapper=text_mapper, property=p2)
            text_actor.position = (0, 0, 0)
            self.scene.add_actor(text_actor)

        else:

            self.scene.add_actor(self.stck_coords[0])
            self.scene.add_actor(self.stck_coords[1])
            self.scene.add_actor(self.stck_coords[2])
            self.scene.add_actor(self.stck_coords[3])
            self.scene.add_actor(self.stck_coords[4])
            self.scene.add_actor(self.stck_coords[5])

            self.scene.mlab.points3d(self.red_balls_data_f1[0],
                                    self.red_balls_data_f1[1],
                                    self.red_balls_data_f1[2], color = (1,0,0),scale_factor=self.red_balls_data_f1[3], resolution=36)


            self.ax = self.scene.mlab.axes(color = (0,0,0), nb_labels = 5, extent = self.axes_limits_f1, line_width=1.0)
            self.ax.axes.font_factor = 0.5


        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)


    def update_plot_vf(self):

        if not self.needUpdate:
            vtext = tvtk.VectorText()
            vtext.text = 'DFormPy 3D'
            text_mapper = tvtk.PolyDataMapper()
            configure_input_data(text_mapper, vtext.get_output())
            vtext.update()
            p2 = tvtk.Property(color=(0, 0.3, 0.3))
            text_actor = tvtk.Follower(mapper=text_mapper, property=p2)
            text_actor.position = (0, 0, 0)
            self.scene.add_actor(text_actor)

        else:

            self.scene.mlab.points3d(self.red_balls_data_vf[0],
                                    self.red_balls_data_vf[1],
                                    self.red_balls_data_vf[2], color = (1,0,0),scale_factor=self.red_balls_data_vf[3], resolution=36)
            
            self.scene.mlab.quiver3d(self.fld_comps[0], self.fld_comps[2], self.fld_comps[4], self.fld_comps[1], self.fld_comps[3], self.fld_comps[5],\
                            colormap='viridis', line_width=3.0, mode='arrow',\
                            scale_factor=0.5, scale_mode = 'vector', opacity = 1)
     

            self.ax = self.scene.mlab.axes(color = (0,0,0), nb_labels = 5, extent = self.axes_limits_vf, line_width=1.0)
            self.ax.axes.font_factor = 0.5

        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)
 
    
    def update_plot_f0(self):

        if not self.needUpdate:
            vtext = tvtk.VectorText()
            vtext.text = 'DFormPy 3D'
            text_mapper = tvtk.PolyDataMapper()
            configure_input_data(text_mapper, vtext.get_output())
            vtext.update()
            p2 = tvtk.Property(color=(0, 0.3, 0.3))
            text_actor = tvtk.Follower(mapper=text_mapper, property=p2)
            text_actor.position = (0, 0, 0)
            self.scene.add_actor(text_actor)

        else:

            cnt = self.scene.mlab.contour3d(self.sc_fld, colormap='jet', opacity = 0.5, contours=7)
            
            self.scene.mlab.colorbar(object = cnt, orientation='vertical')

            self.scene.mlab.outline(line_width=1.0)
            self.ax = self.scene.mlab.axes(color = (0,0,0), ranges=(lft,rght,lft,rght,lft,rght), nb_labels = 5, line_width=1.0)
            self.ax.axes.font_factor = 0.5

        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)
 

    def update_plot_f2(self):

        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)

        if not self.needUpdate:
            vtext = tvtk.VectorText()
            vtext.text = 'DFormPy 3D'
            text_mapper = tvtk.PolyDataMapper()
            configure_input_data(text_mapper, vtext.get_output())
            vtext.update()
            p2 = tvtk.Property(color=(0, 0.3, 0.3))
            text_actor = tvtk.Follower(mapper=text_mapper, property=p2)
            text_actor.position = (0, 0, 0)
            self.scene.add_actor(text_actor)

        else:

            f2 = df3.form_2_3d(self.x, self.y, self.z, self.fx, self.fy, self.fz)
            f2.give_eqn(self.fx_eqn, self.fy_eqn, self.fz_eqn)
            f2.plot(self.scene)

 
    def update_plot_f3(self):

        self.scene.background = (1, 1, 1)
        self.scene.foreground = (0, 0, 0)

        if not self.needUpdate:
            vtext = tvtk.VectorText()
            vtext.text = 'DFormPy 3D'
            text_mapper = tvtk.PolyDataMapper()
            configure_input_data(text_mapper, vtext.get_output())
            vtext.update()
            p2 = tvtk.Property(color=(0, 0.3, 0.3))
            text_actor = tvtk.Follower(mapper=text_mapper, property=p2)
            text_actor.position = (0, 0, 0)
            self.scene.add_actor(text_actor)

        else:

            f3 = df3.form_3_3d(self.x3, self.y3, self.z3, self.f3)
            f3.give_eqn(self.f3_eqn)
            f3.plot(self.scene)



 







 
################################################################################

 
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        
        layout = QtGui.QGridLayout(self)

        layout.setSpacing(0)
        self.visualization = Visualization()

        spc = QtGui.QLabel('')


        self.label1 = QtGui.QLabel('dx')
        self.label1.setAlignment(Qt.AlignCenter)
        self.line_edit1 = QtGui.QLineEdit("x/sqrt(x**2+y**2-z**2)")
        self.box1 = QCheckBox(" ")
        self.box1.setEnabled(False)

        self.label2 = QtGui.QLabel('dy')
        self.label2.setAlignment(Qt.AlignCenter)
        self.line_edit2 = QtGui.QLineEdit("y/sqrt(x**2+y**2-z**2)")
        self.box2 = QCheckBox(" ")
        self.box2.setEnabled(False)

        self.label3 = QtGui.QLabel('dz')
        self.label3.setAlignment(Qt.AlignCenter)
        self.line_edit3 = QtGui.QLineEdit("z/sqrt(x**2+y**2-z**2)")
        self.box3 = QCheckBox(" ")
        self.box3.setEnabled(False)

        self.combobox1 = QtGui.QComboBox()
        self.combobox1.addItems(['Covariant Field (Vector Field)', '0-form (Scalar Field)', '1-form', '2-form', '3-form'])
        self.combobox1.setCurrentIndex(2)



        self.sublayout1 = QtGui.QHBoxLayout()
        self.sublayout1.addWidget(self.line_edit1, 9)
        self.sublayout1.addWidget(self.box1, 1, alignment=Qt.AlignRight)

        self.sublayout2 = QtGui.QHBoxLayout()
        self.sublayout2.addWidget(self.line_edit2, 9)
        self.sublayout2.addWidget(self.box2, 1, alignment=Qt.AlignRight)

        self.sublayout3 = QtGui.QHBoxLayout()
        self.sublayout3.addWidget(self.line_edit3, 9)
        self.sublayout3.addWidget(self.box3, 1, alignment=Qt.AlignRight)
        
        

        

        def disableWidget():
            if self.combobox1.currentIndex() == 0:
                self.label1.setText('x component')
                self.label2.setText('y component')
                self.label3.setText('z component')
                self.line_edit1.setEnabled(True)
                self.line_edit2.setEnabled(True)
                self.line_edit3.setEnabled(True)
                self.label1.setEnabled(True)
                self.label2.setEnabled(True)
                self.label3.setEnabled(True)
                self.box1.setEnabled(False)
                self.box2.setEnabled(False)
                self.box3.setEnabled(False)

            if self.combobox1.currentIndex() == 1:
                self.label1.setText('Field')
                self.label2.setText(' ')
                self.label3.setText(' ')
                self.line_edit1.setEnabled(True)
                self.line_edit2.setEnabled(False)
                self.line_edit3.setEnabled(False)
                self.label1.setEnabled(True)
                self.label2.setEnabled(False)
                self.label3.setEnabled(False)
                self.box1.setEnabled(False)
                self.box2.setEnabled(False)
                self.box3.setEnabled(False)

            
            if self.combobox1.currentIndex() == 2:
                self.label1.setText('dx')
                self.label2.setText('dy')
                self.label3.setText('dz')
                self.line_edit1.setEnabled(True)
                self.line_edit2.setEnabled(True)
                self.line_edit3.setEnabled(True)
                self.label1.setEnabled(True)
                self.label2.setEnabled(True)
                self.label3.setEnabled(True)
                self.box1.setEnabled(False)
                self.box2.setEnabled(False)
                self.box3.setEnabled(False)

            if self.combobox1.currentIndex() == 3:
                self.label1.setText('dy ∧ dz')
                self.label2.setText('dz ∧ dx')
                self.label3.setText('dx ∧ dy')
                self.line_edit1.setEnabled(False)
                self.line_edit2.setEnabled(False)
                self.line_edit3.setEnabled(False)
                self.label1.setEnabled(True)
                self.label2.setEnabled(True)
                self.label3.setEnabled(True)
                self.box1.setEnabled(True)
                self.box2.setEnabled(True)
                self.box3.setEnabled(True)

            if self.combobox1.currentIndex() == 4:
                self.label1.setText('dx ∧ dy ∧ dz')
                self.label2.setText(' ')
                self.label3.setText(' ')
                self.line_edit1.setEnabled(True)
                self.line_edit2.setEnabled(False)
                self.line_edit3.setEnabled(False)
                self.label1.setEnabled(True)
                self.label2.setEnabled(False)
                self.label3.setEnabled(False)
                self.box1.setEnabled(False)
                self.box2.setEnabled(False)
                self.box3.setEnabled(False)







        self.combobox1.currentIndexChanged['QString'].connect(disableWidget)

        self.box2.toggled.connect(self.box1.setDisabled)
        self.box2.toggled.connect(self.box3.setDisabled)
        self.box2.toggled.connect(self.line_edit2.setEnabled)

        self.box1.toggled.connect(self.box2.setDisabled)
        self.box1.toggled.connect(self.box3.setDisabled)
        self.box1.toggled.connect(self.line_edit1.setEnabled)

        self.box3.toggled.connect(self.box1.setDisabled)
        self.box3.toggled.connect(self.box2.setDisabled)
        self.box3.toggled.connect(self.line_edit3.setEnabled)



        self.ui = self.visualization.edit_traits(parent=self, 
                                                 kind='subpanel').control
    
        layout.addWidget(self.ui, 0, 0)
        layout.addWidget(spc, 1, 0)
        layout.addWidget(self.combobox1, 2, 0)


        layout.addWidget(self.label1, 3, 0)
        layout.addLayout(self.sublayout1, 4, 0)


        layout.addWidget(self.label2, 5, 0)
        layout.addLayout(self.sublayout2, 6, 0)


        layout.addWidget(self.label3, 7, 0)
        layout.addLayout(self.sublayout3, 8, 0)



        

        self.ui.setParent(self)
 

    def clear(self):
        self.visualization.clear()




    def create_df3_plot(self):

        if self.combobox1.currentIndex()==2:
            stck_coords, red_balls_data, axes_limits = self.createForm1()
            self.visualization.takePlotParametresF1(stck_coords, red_balls_data, axes_limits)
        elif self.combobox1.currentIndex()==0:
            fld_copms, red_balls_data1, axes_limits1 = self.createVF()
            self.visualization.takePlotParametresVF(fld_copms, red_balls_data1, axes_limits1)
        elif self.combobox1.currentIndex()==1:
            sc_fld, axes_limits_f0 = self.createForm0()
            self.visualization.takePlotParametresF0(sc_fld, axes_limits_f0)
        elif self.combobox1.currentIndex()==3:
            xg, yg, zg, fx, fy, fz, fx_eqn, fy_eqn, fz_eqn = self.createForm2()
            self.visualization.takePlotParametresF2(xg, yg, zg, fx, fy, fz, fx_eqn, fy_eqn, fz_eqn)
        elif self.combobox1.currentIndex()==4:
            xg3, yg3, zg3, f3, f3_eqn = self.createForm3()
            self.visualization.takePlotParametresF3(xg3, yg3, zg3, f3, f3_eqn)
        else:
            print('kek')
 


    def createForm1(self):

        grid = np.linspace(lft, rght, pts)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        fx = xg/np.sqrt(xg**2+yg**2 - zg**2)
        fy = yg/np.sqrt(xg**2+yg**2- zg**2)
        fz = yg/np.sqrt(xg**2+yg**2- zg**2)
        

        fx_eqn = self.line_edit1.text()
        fy_eqn = self.line_edit2.text()
        fz_eqn = self.line_edit3.text()

        form_1 = df3.form_1_3d(xg, yg, zg, fx, fy, fz)

        form_1.give_eqn(fx_eqn, fy_eqn, fz_eqn)

        stck_coords, red_balls_data, axes_limits = form_1.plot()



        return stck_coords, red_balls_data, axes_limits
 

    def createForm0(self):

        grid = np.linspace(lft, rght, pts)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        f0 = xg/np.sqrt(xg**2+yg**2 - zg**2)

        f0_eqn = self.line_edit1.text()

        form_0 = df3.form_0_3d(xg, yg, zg, f0)

        form_0.give_eqn(f0_eqn)

        sc_field, axes_limits = form_0.plot()



        return sc_field, axes_limits


    def createVF(self):

        grid = np.linspace(lft, rght, pts)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        fx = xg/np.sqrt(xg**2+yg**2 - zg**2)
        fy = yg/np.sqrt(xg**2+yg**2- zg**2)
        fz = yg/np.sqrt(xg**2+yg**2- zg**2)
        

        fx_eqn = self.line_edit1.text()
        fy_eqn = self.line_edit2.text()
        fz_eqn = self.line_edit3.text()

        vf = df3.vector_field3(xg, yg, zg, fx, fy, fz)

        vf.give_eqn(fx_eqn, fy_eqn, fz_eqn)

        fld_comps, red_balls_data, axes_limits = vf.plot()



        return fld_comps, red_balls_data, axes_limits

    
    def createForm2(self):

        grid = np.linspace(lft, rght, pts)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        fx = xg/np.sqrt(xg**2+yg**2 - zg**2)
        fy = yg/np.sqrt(xg**2+yg**2- zg**2)
        fz = yg/np.sqrt(xg**2+yg**2- zg**2)
        

        

        if self.box1.checkState()==2:
            fx_eqn = self.line_edit1.text()
            fy_eqn = '0'
            fz_eqn = '0'

        if self.box2.checkState()==2:
            fx_eqn = '0'
            fy_eqn = self.line_edit2.text()
            fz_eqn = '0'

        if self.box3.checkState()==2:
            fx_eqn = '0'
            fy_eqn = '0'
            fz_eqn = self.line_edit3.text()

        if self.box1.checkState()==0 and self.box2.checkState()==0 and self.box3.checkState()==0:
            fx_eqn = '1'
            fy_eqn = '0'
            fz_eqn = '0'

            self.line_edit1.setText('1')
            self.line_edit2.setText('0')
            self.line_edit3.setText('0')



        return xg, yg, zg, fx, fy, fz, fx_eqn, fy_eqn, fz_eqn


    
    def createForm3(self):

        grid = np.linspace(lft, rght, pts)

        xg, yg, zg = np.meshgrid(grid, grid, grid)

        f3 = xg/np.sqrt(xg**2+yg**2 - zg**2)
        
        f3_eqn = self.line_edit1.text()

        return xg, yg, zg, f3, f3_eqn







if __name__ == "__main__":

    app = QtGui.QApplication.instance()
 
    container = QtGui.QWidget()
    container.setWindowTitle("DFormPy 3D GUI")

    overlayout = QtGui.QHBoxLayout(container)
    container.setGeometry(0,0,500,500)

    layout = QtGui.QVBoxLayout()
    
 
    mayavi_widget = MayaviQWidget()
 
    button = QtGui.QPushButton('Generate plot')
    button.clicked.connect(mayavi_widget.create_df3_plot)

    
 
    layout.addWidget(mayavi_widget)
    
    layout.addWidget(button)



    button2 = QtGui.QPushButton('Clear canvas')
    button2.move(100,0)
    button2.clicked.connect(mayavi_widget.clear)
    overlayout.addLayout(layout)
    overlayout.addWidget(button2)
    
 
    container.show()
 
    app.exec_()



