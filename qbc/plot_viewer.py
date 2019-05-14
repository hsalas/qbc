'''Creates GUI for ploting the results of QbC_mgii_v6
'''
import sys
import os
 
from xkcd_rgb import *
 
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
 
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from astropy.table import Table, QTable, Column
from plots import plot_dNdz, plot_dNdz_vs_x
import matplotlib.pyplot as plt
import random
import pickle
 
# import resources
# import resources2
# import resources_local_b12
# import resources_local_b40
 
# matplotlib style
plt.rcParams['figure.figsize'] = [4.0, 3.0]
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plt.rcParams['xtick.major.size'] = 12.0
plt.rcParams['ytick.major.size'] = 12.0
plt.rcParams['xtick.minor.size'] = 6.0
plt.rcParams['ytick.minor.size'] = 6.0
 
 
#parmaters to use
 
b = [12, 40]
n = [0.5, 1.0, 2.0] 
sig = [0.0, 1.0, 2.0, 3.0]
ewlims = [(0.6, 2.0), (0.6, 1.0), (1.0, 1.5), (1.5, 2.0) ]
masslims = [ (13.6, 16.0), (13.6, 14.0), (14.0, 14.2), (14.20, 16.0)]
zlims = [('0.36','0.44'), ('0.44','0.52'), ('0.52','0.60')]
distances = ['com', 'pro', 'r200']
grids = ['log', 'snp', 'linear']
cl_types = ['spec','phot']
signal2noise = ['local', '']
limit_by = ['lya', 'civ']
color_counter = 0
plots = []
 
class MatplolibPlot(QDialog):
    """Creates a Widget containing a plot using matplotlib"""
    def __init__(self):
        super(MatplolibPlot, self).__init__()
 
        # a figure and ax instance to plot on
        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
         
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
 
        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        # self.toolbar.setFixedHeight(20)
 
        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        # self.button.clicked.connect(self.plot)
 
        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)
 
        # appearance
        self.labels_size = 22
        self.legend_size =16
        self.annotate_size = 12
        self.ms = 14
        self.marker = 'o'
        self.xlim = (0.1, 40)
        self.ylim = (0.01, 30)
        self.elinewidth = 2
 
    def make_plot(self, X, table, color, legend):
        # plot_dNdz(self.ax, table , marker=self.marker, annotate=True, c=color, capsize=0, lw=0.5, ms=self.ms, ls='None', elinewidth=self.elinewidth, annotate_size=self.annotate_size)
        plot_dNdz_vs_x(X, self.ax, table , marker=self.marker, annotate=True, c=color, capsize=0, lw=0.5, ms=self.ms, ls='None', elinewidth=self.elinewidth, annotate_size=self.annotate_size)
        self.ax.set_ylabel('dn/dz', fontsize=self.labels_size)
        self.ax.set_xlabel('Impact parameter', fontsize=self.labels_size)
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.legend(legend, loc='upper right', fontsize=self.legend_size, numpoints=1)   
        self.ax.minorticks_on()
        self.canvas.draw()
 
    def plot_field_value(self, field_value):
        self.ax.plot(self.xlim, [field_value]*2, 'k--', lw=2, label='_nolegend_')
        self.canvas.draw()
 
    def clear_plots(self):
        self.figure.clf()
        self.ax = self.figure.add_subplot(111)
        self.canvas.draw()
        del plots[:]
 
class RadioButtonWidget(QWidget):
    '''this class cretes a gropu of radio buttons froam a given list of labels'''
 
    #constructor
    def __init__(self, label, instruction, button_list):
        super(RadioButtonWidget, self).__init__()#call the super class constructor
 
        #create widgets
        self.title_label = QLabel(label)
        self.radio_group_box = QGroupBox(instruction)
        self.radio_button_group = QButtonGroup()
 
        #set horizontal alingment for label
        self.title_label.setAlignment(Qt.AlignHCenter)
 
        #create the radio buttons
        self.radio_button_list = []
        for each in button_list:
            each_s = str(each)
            self.radio_button_list.append(QRadioButton(each_s))
 
        #set the default checked value
        self.radio_button_list[0].setChecked(True)
 
        #create layout for radio buttons and add them
        self.radio_button_layout = QVBoxLayout()
 
        #add buttons to the layout  and button group
        counter = 0
        for each in self.radio_button_list:
            self.radio_button_layout.addWidget(each)
            self.radio_button_group.addButton(each)
            self.radio_button_group.setId(each, counter)
            counter +=1
 
        #add radio buttons to the group
        self.radio_group_box.setLayout(self.radio_button_layout)
 
        #create a layout for the whole widget
        self.main_layout = QVBoxLayout()
 
        #add widgets to the layout
        self.main_layout.addWidget(self.title_label)
        self.main_layout.addWidget(self.radio_group_box)
 
        #set aligment for the layout
        self.main_layout.setAlignment(Qt.AlignTop)
 
        #set the layout for this widget
        self.setLayout(self.main_layout)
 
    #method to find out the selected button
    def selected_button(self):
        return self.radio_button_group.checkedId()
 
    #method for change default vaue
    def select_default(self, index):
        self.radio_button_list[index].setChecked(True)
         
class LineEdit(QWidget):
    """Creates a line edit widget that accepts floats, given a title and an instruction"""
    def __init__(self, title, instruction):
        super(LineEdit, self).__init__()
        self.title_label = QLabel(title)
        self.lineedit_groupbox = QGroupBox(instruction)
        self.lineedit = QLineEdit()
 
        #set up a validator to ensure only real numbers can be entered
        self.validator = QDoubleValidator()
        self.lineedit.setValidator(self.validator)
 
        #create the layout for the widget
        self.lineedit_layout = QHBoxLayout()
 
        #add lineedit to groupbox
        self.lineedit_layout.addWidget(self.lineedit)
        self.lineedit_groupbox.setLayout(self.lineedit_layout)
 
        #create a layout for the whole widget
        self.main_layout = QVBoxLayout()
 
        #add widgets to the layout
        self.main_layout.addWidget(self.title_label)
        self.main_layout.addWidget(self.lineedit_groupbox)
 
        #set aligment for the layout
        self.main_layout.setAlignment(Qt.AlignTop)
 
        #set the layout for this widget
        self.setLayout(self.main_layout)
 
    def get_number(self):
        num = self.lineedit.text()
        num = float(num)
        return (num)
     
    def set_default(self, num):
        text = '{:.1f}'.format(num)
        self.lineedit.setText(text)
 
class SpinBox(QWidget):
    """Creates a spinbox from a label, a maximum and a minumum"""
 
    def __init__(self, label, min_value, max_value, suffix='', default_value=0):
        super(SpinBox, self).__init__()#call the super class constructor
 
        #create widgets
        self.title_label = QLabel(label)
        self.space_label =QLabel('')
        self.spinbox = QSpinBox()
        self.spinbox.setRange(min_value, max_value)
        self.spinbox.setSuffix(suffix)
        self.spinbox.setValue(default_value)
 
        #adjust size and alingment for label widget
        # self.title_label.setFixedHeight(45)
        self.spinbox.setFixedWidth(100)
         
        #create a layout for the whole widgrt
        self.main_layout = QVBoxLayout()
 
        #add widgets to layout
        self.main_layout.addWidget(self.title_label)
        self.main_layout.addWidget(self.space_label)
        self.main_layout.addWidget(self.spinbox)
 
        #set aligment for the layout
        self.main_layout.setAlignment(Qt.AlignTop)
 
        #set the layout for this widget
        self.setLayout(self.main_layout)
 
    #method to find out the selected button
    def value(self):
        return int(self.spinbox.value())
 
class ListWidget(QWidget):
    """Creates a list widget from a labela, a list of options pairs (a, b), and a units"""
    def __init__(self, label, options_list, units):
        super(ListWidget, self).__init__()
 
        #create widgets
        self.title_label = QLabel(label+' '+units)
        self.listwidget = QListWidget()
 
        #add options to the list widget
        for each in options_list:
            s_min = '{}'.format(each[0])
            s_max = '{}'.format(each[1])
            s = s_min+' to '+s_max#+' '+units
            self.item =QListWidgetItem(s)
            self.listwidget.addItem(self.item)
 
        #select default value
        self.listwidget.setCurrentRow(0)
 
        #adjust size and aligment for widgets
        self.title_label.setFixedHeight(35)
        self.title_label.setAlignment(Qt.AlignHCenter)
        self.title_label.setAlignment(Qt.AlignTop)
        self.listwidget.setFixedSize(self.listwidget.sizeHintForColumn(0) + 2 * self.listwidget.frameWidth(), self.listwidget.sizeHintForRow(0) * self.listwidget.count() + 2 * self.listwidget.frameWidth())
         
        #create layout for the whole widget
        self.main_layout = QVBoxLayout()
 
        #add widgets to the layout
        self.main_layout.addWidget(self.title_label)
        self.main_layout.addWidget(self.listwidget)
 
        #set aligment for the layout
        self.main_layout.setAlignment(Qt.AlignTop)
 
        #set the layout for this widget
        self.setLayout(self.main_layout)
 
    #metod to find out selected option
    def selected(self):
        return self.listwidget.currentRow(), self.listwidget.currentItem().text()
         
class PlotView(QGraphicsView):
     
    #constructor
    def __init__(self):
        super(PlotView, self).__init__()
         
        #create view
        self.view = QGraphicsView() 
        self.scene = QGraphicsScene(self.view)
        self.view.setAlignment(Qt.AlignCenter)
        self.scene.update()
 
        self.label = QLabel()
        self.label.setAlignment(Qt.AlignCenter)
         
        #create a layout for the whole widgrt
        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.view)
 
    def load_plot(self, name):
        #get the graphics
        self.plot = QPixmap(':/{0}'.format(name))
 
        #add the graphics to scene 
        self.plot = self.plot.scaled(1500, 750, Qt.KeepAspectRatio)
         
        self.rect = QRect(150, 40, 1210 , 680);
        self.cropped = self.plot.copy(self.rect)        
         
        self.scene.addPixmap(self.cropped)
        self.setScene(self.scene)
         
        #update plot
        self.scene.update()
 
 
class PlotWindow(QMainWindow):
    '''this creates a main window to create the plots from QbC_mgii_v7.py using matplotlib'''
 
    #constructor
    def __init__(self):
        super(PlotWindow, self).__init__()#call super class constructor
        self.create_select_plot_layout()
 
    def create_select_plot_layout(self):
        #this is the initial layout of the window - to select the plot parameter
        self.height_select_area = 120
        self.width_select_area = 200
        self.field_value = 0.3
        self.cleared = 'no'
        self.distance_default_value = 0
 
        #radio buttons for impact parameter, distance, grid, cluster type, limit
        self.impact_parameter_buttons =RadioButtonWidget('Impact Parameter', '', b)
        self.distance_buttons = RadioButtonWidget('Distance Units', '', ['Mpc Comoving', 'Mpc Proper', 'R 200'])
        self.cluster_type_buttons = RadioButtonWidget('Cluster Z Type', '', ['Spectroscopic','Photometric'])
        self.grid_buttons = RadioButtonWidget('Grid', '', ['Log', 'SNP', 'Linear'])
        self.signal2noise_buttons = RadioButtonWidget('S2N', '', ['Local', 'Global'])
        self.limit_by_buttons = RadioButtonWidget('Limit by', '', ['Lyman alpha', 'CIV'])
        self.nbins_widget = RadioButtonWidget('#bins', '', ['x0.5', 'x1.0', 'x2.0'])
        self.impact_parameter_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.distance_buttons.setFixedSize(self.width_select_area, self.height_select_area + 15)
        self.grid_buttons.setFixedSize(self.width_select_area / 2.0, self.height_select_area + 15)
        self.nbins_widget.setFixedSize(self.width_select_area /2.0, self.height_select_area + 15)
        self.cluster_type_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.signal2noise_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.limit_by_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.nbins_widget.select_default(1)
 
        #signifficance spinbox, n spinboc
        self.significance_spinbox = SpinBox('Significance', 0, 3)
        self.significance_spinbox.setFixedSize(self.width_select_area, self.height_select_area)
 
        #list widgets for Mass and EW
        self.mass_listwiget = ListWidget('Mass', masslims, '[log(Msun)]')
        self.ew_listwidget = ListWidget('EW', ewlims, '[AA]')
        self.mass_listwiget.setFixedSize(self.width_select_area, self.height_select_area)
        self.ew_listwidget.setFixedSize(self.width_select_area, self.height_select_area)
 
        #create field input box
        self.field_widget = LineEdit('Field Value', 'Insert field value')
        self.field_widget.setFixedSize(self.width_select_area, self.height_select_area)
        self.field_widget.set_default(self.field_value)
 
        #Plot and Clear Buttons 
        self.plot_button = QPushButton('Plot')
        self.plot_button.setFixedSize(self.width_select_area, 30)
        self.clear_button = QPushButton('Clear')
        self.clear_button.setFixedSize(self.width_select_area, 30)
 
        #create plot from matplotlib
        self.plot_view = MatplolibPlot()
        # self.plot_view.sizeHint=500       
 
        #create layouts
        self.grid_layout = QHBoxLayout()
        self.options_layout = QGridLayout()
        self.graphic_layout = QHBoxLayout()
        self.buttons_layout = QVBoxLayout()
 
        #add widgets to grid layput
        self.grid_layout.addWidget(self.grid_buttons)
        self.grid_layout.addWidget(self.nbins_widget)
 
        #add widgets/layput to the options layout
        self.options_layout.addWidget(self.impact_parameter_buttons, 0, 0)
        self.options_layout.addWidget(self.cluster_type_buttons, 0,1)
        self.options_layout.addWidget(self.signal2noise_buttons,1, 0 )
        self.options_layout.addWidget(self.limit_by_buttons, 1, 1)
        self.options_layout.addLayout(self.grid_layout, 2, 0)
        self.options_layout.addWidget(self.distance_buttons, 2, 1)
        self.options_layout.addWidget(self.mass_listwiget, 3, 0)
        self.options_layout.addWidget(self.ew_listwidget, 3, 1)
        self.options_layout.addWidget(self.field_widget, 4, 0)
        self.options_layout.addWidget(self.significance_spinbox, 4, 1)
        self.options_layout.addWidget(self.plot_button, 5, 0)
        self.options_layout.addWidget(self.clear_button, 5,1)
 
        #add widget/layout to the graphics layout
        self.graphic_layout.addLayout(self.options_layout)
        self.graphic_layout.addWidget(self.plot_view)
 
        #create plot widget and add it to layout
        self.select_plot_widget = QWidget()
        self.select_plot_widget.setLayout(self.graphic_layout)
 
        self.setCentralWidget(self.select_plot_widget)
 
        #connections
        self.plot_button.clicked.connect(self.instantiate_plot)
        self.clear_button.clicked.connect(self.clear_plot)
 
    def update_plot(self, X, name, color, legend):
        #update the plot image
        self.plot_view.make_plot(X, name, color, legend)
 
    def plot_field(self, field):
        #plots the field value
        self.plot_view.plot_field_value(field)
 
    def clear_plot(self):
        self.plot_view.clear_plots()
        self.cleared = 'yes'
 
    def instantiate_plot(self):
 
        impact_parameter_value = self.impact_parameter_buttons.selected_button() #get the radio that was selected
        impact_parameter = b[impact_parameter_value]    
        grid_value = self.grid_buttons.selected_button()
        grid = grids[grid_value]    
        # print('Grid: {}'.format(grid))
 
        cluster_type_value = self.cluster_type_buttons.selected_button()
        cl_type = cl_types[cluster_type_value]
        # print('Cluster  type: {}'.format(cl_type))
         
        distance_value = self.distance_buttons.selected_button()
        distance = distances[distance_value]
        # print('Distance units: {}'.format(distance))
        if distance_value != self.distance_default_value:
            self.clear_plot()
            self.distance_default_value = distance_value
 
        signal2noise_value = self.signal2noise_buttons.selected_button()
        s2n = signal2noise[signal2noise_value]
 
        significance_value = self.significance_spinbox.value()
        # print('Significance: {}'.format(significance_value))
 
        mass_lims_index, mass_lims = self.mass_listwiget.selected()
        masslim = masslims[mass_lims_index]
        # print('{}, index: {}'.format(mass_lims, mass_lims_index))
 
        ew_lims_index, ew_lims = self.ew_listwidget.selected()
        ewlim = ewlims[ew_lims_index]
        # print('{}, index: {}'.format(ew_lims, ew_lims_index))
 
        limit_by_value = self.limit_by_buttons.selected_button()
        limit = limit_by[limit_by_value]
 
        nbins_value = self.nbins_widget.selected_button()
        nbins = n[nbins_value]
 
        field =self.field_widget.get_number()
        # print('{:.1f}'.format(field))
 
        zlim=(0,10)
 
        grid_str = '{}_b{:.1f}_n{:.1f}'.format(grid, impact_parameter, nbins)
        alias = 'mass {} to {}, ew {} to {}, s {}, {}, {}'.format(masslim[0], masslim[1], ewlim[0], ewlim[1], significance_value, cl_type, limit)
        # dir_name = '../saved_files/{}-lim_{}-mass_10e{}_to_10e{}-rew_{}_to_{}-s_{:.1f}_{}/{}/results_{}.pickle'.format(grid_str,limit, masslim[0], masslim[1], ewlim[0], ewlim[1] , significance_value ,s2n, cl_type, distance)
        dir_name = '../saved_files/dndz_v_b/{}-lim_{}-mass_10e{}_to_10e{}-rew_{}_to_{}-z_{:.2f}_to_{:.2f}-s_{:.1f}_{}/{}/results_{}.pickle'.format(grid_str,limit, masslim[0], masslim[1], ewlim[0], ewlim[1], zlim[0], zlim[1], significance_value ,s2n, cl_type, distance)
        color_list = xkcd_color_list(banned=['white'])
        random.shuffle(color_list)
        # print(dir_name)
         
        if os.path.isfile(dir_name):
            try:
                with open(dir_name, 'rb') as f:
                    results = pickle.load(f)
            except UnicodeDecodeError:
                with open(dir_name, 'rb') as f:
                    results = pickle.load(f, encoding='latin1')
            plots.append(alias)
            self.update_plot('b', results, color_list[len(plots) - 1], plots)
            if len(plots) ==1:
                self.plot_field(field)
            elif field != self.field_value:
                self.plot_field(field)
                self.field_value = field
            elif self.cleared == 'yes':
                self.plot_field(field)
                self.cleared = 'no'
        else:
            print('{} File not found, make sure QbC_mgii_v6.py was run with this configuration'.format(dir_name))
             
 
class PlotWindowZ(QMainWindow):
    '''this creates a main window to create the plots from QbC_mgii_v6.py using matplotlib'''
 
    #constructor
    def __init__(self):
        super(PlotWindowZ, self).__init__()#call super class constructor
        self.create_select_plot_layout()
 
    def create_select_plot_layout(self):
        #this is the initial layout of the window - to select the plot parameter
        self.height_select_area = 120
        self.width_select_area = 200
        self.field_value = 0.3
        self.cleared = 'no'
        self.distance_default_value = 0
 
        #radio buttons for impact parameter, distance, grid, cluster type, limit
        self.impact_parameter_buttons =RadioButtonWidget('Impact Parameter', '', b)
        self.distance_buttons = RadioButtonWidget('Distance Units', '', ['Mpc Comoving', 'Mpc Proper', 'R 200'])
        self.cluster_type_buttons = RadioButtonWidget('Cluster Z Type', '', ['Spectroscopic','Photometric'])
        self.grid_buttons = RadioButtonWidget('Grid', '', ['Log', 'SNP', 'Linear'])
        self.signal2noise_buttons = RadioButtonWidget('S2N', '', ['Local', 'Global'])
        self.limit_by_buttons = RadioButtonWidget('Limit by', '', ['Lyman alpha', 'CIV'])
        self.nbins_widget = RadioButtonWidget('#bins', '', ['x0.5', 'x1.0', 'x2.0'])
        self.impact_parameter_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.distance_buttons.setFixedSize(self.width_select_area, self.height_select_area + 15)
        self.grid_buttons.setFixedSize(self.width_select_area / 2.0, self.height_select_area + 15)
        self.nbins_widget.setFixedSize(self.width_select_area /2.0, self.height_select_area + 15)
        self.cluster_type_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.signal2noise_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.limit_by_buttons.setFixedSize(self.width_select_area, self.height_select_area)
        self.nbins_widget.select_default(1)
 
        #signifficance spinbox, n spinboc
        self.significance_spinbox = SpinBox('Significance', 0, 3)
        self.significance_spinbox.setFixedSize(self.width_select_area, self.height_select_area)
 
        #list widgets for Mass and EW
        self.mass_listwiget = ListWidget('Mass', masslims, '[log(Msun)]')
        self.ew_listwidget = ListWidget('EW', ewlims, '[AA]')
        self.z_listwidget = ListWidget('Z', zlims, '')
        self.mass_listwiget.setFixedSize(self.width_select_area, self.height_select_area)
        self.ew_listwidget.setFixedSize(self.width_select_area, self.height_select_area)
        self.z_listwidget.setFixedSize(self.width_select_area, self.height_select_area)
 
        #create field input box
        self.field_widget = LineEdit('Field Value', 'Insert field value')
        self.field_widget.setFixedSize(self.width_select_area, self.height_select_area)
        self.field_widget.set_default(self.field_value)
 
        #Plot and Clear Buttons 
        self.plot_button = QPushButton('Plot')
        self.plot_button.setFixedSize(self.width_select_area, 30)
        self.clear_button = QPushButton('Clear')
        self.clear_button.setFixedSize(self.width_select_area, 30)
 
        #create plot from matplotlib
        self.plot_view = MatplolibPlot()
        # self.plot_view.sizeHint=500       
 
        #create layouts
        self.grid_layout = QHBoxLayout()
        self.options_layout = QGridLayout()
        self.graphic_layout = QHBoxLayout()
        self.buttons_layout = QVBoxLayout()
 
        #add widgets to grid layput
        self.grid_layout.addWidget(self.grid_buttons)
        self.grid_layout.addWidget(self.nbins_widget)
 
        #add widgets/layput to the options layout
        self.options_layout.addWidget(self.impact_parameter_buttons, 0, 0)
        self.options_layout.addWidget(self.cluster_type_buttons, 0,1)
        self.options_layout.addWidget(self.z_listwidget ,1, 0 )
        self.options_layout.addWidget(self.limit_by_buttons, 1, 1)
        self.options_layout.addLayout(self.grid_layout, 2, 0)
        self.options_layout.addWidget(self.distance_buttons, 2, 1)
        self.options_layout.addWidget(self.mass_listwiget, 3, 0)
        self.options_layout.addWidget(self.ew_listwidget, 3, 1)
        self.options_layout.addWidget(self.field_widget, 4, 0)
        self.options_layout.addWidget(self.significance_spinbox, 4, 1)
        self.options_layout.addWidget(self.plot_button, 5, 0)
        self.options_layout.addWidget(self.clear_button, 5,1)
 
        #add widget/layout to the graphics layout
        self.graphic_layout.addLayout(self.options_layout)
        self.graphic_layout.addWidget(self.plot_view)
 
        #create plot widget and add it to layout
        self.select_plot_widget = QWidget()
        self.select_plot_widget.setLayout(self.graphic_layout)
 
        self.setCentralWidget(self.select_plot_widget)
 
        #connections
        self.plot_button.clicked.connect(self.instantiate_plot)
        self.clear_button.clicked.connect(self.clear_plot)
 
    def update_plot(self, name, color, legend):
        #update the plot image
        self.plot_view.make_plot(name, color, legend)
 
    def plot_field(self, field):
        #plots the field value
        self.plot_view.plot_field_value(field)
 
    def clear_plot(self):
        self.plot_view.clear_plots()
        self.cleared = 'yes'
 
    def instantiate_plot(self):
 
        impact_parameter_value = self.impact_parameter_buttons.selected_button() #get the radio that was selected
        impact_parameter = b[impact_parameter_value]    
        grid_value = self.grid_buttons.selected_button()
        grid = grids[grid_value]    
        # print('Grid: {}'.format(grid))
 
        cluster_type_value = self.cluster_type_buttons.selected_button()
        cl_type = cl_types[cluster_type_value]
        # print('Cluster  type: {}'.format(cl_type))
         
        distance_value = self.distance_buttons.selected_button()
        distance = distances[distance_value]
        # print('Distance units: {}'.format(distance))
        if distance_value != self.distance_default_value:
            self.clear_plot()
            self.distance_default_value = distance_value
 
        signal2noise_value = self.signal2noise_buttons.selected_button()
        s2n = signal2noise[signal2noise_value]
 
        significance_value = self.significance_spinbox.value()
        # print('Significance: {}'.format(significance_value))
 
        mass_lims_index, mass_lims = self.mass_listwiget.selected()
        masslim = masslims[mass_lims_index]
        # print('{}, index: {}'.format(mass_lims, mass_lims_index))
 
        ew_lims_index, ew_lims = self.ew_listwidget.selected()
        ewlim = ewlims[ew_lims_index]
        # print('{}, index: {}'.format(ew_lims, ew_lims_index))
 
        limit_by_value = self.limit_by_buttons.selected_button()
        limit = limit_by[limit_by_value]
 
        nbins_value = self.nbins_widget.selected_button()
        nbins = n[nbins_value]
 
        zlims_index, z_lims = self.z_listwidget.selected()
        zlim = zlims[zlims_index]
 
        field =self.field_widget.get_number()
        # print('{:.1f}'.format(field))
 
        grid_str = '{}_b{:.1f}_n{:.1f}'.format(grid, impact_parameter, nbins)
        alias = 'mass {} to {}, ew {} to {}, z {} to {}, s {}, {}, {}'.format(masslim[0], masslim[1], ewlim[0], ewlim[1], zlim[0], zlim[1], significance_value, cl_type, limit)
        dir_name = '../saved_files/{}-lim_{}-mass_10e{}_to_10e{}-rew_{}_to_{}-z_{}_to_{}-s_{:.1f}_{}/{}/results_{}.pickle'.format(grid_str,limit, masslim[0], masslim[1], ewlim[0], ewlim[1], zlim[0], zlim[1], significance_value ,s2n, cl_type, distance)
        color_list = xkcd_color_list(banned=['white'])
        random.shuffle(color_list)
        # print(dir_name)
         
        if os.path.isfile(dir_name):
            try:
                with open(dir_name, 'rb') as f:
                    results = pickle.load(f)
            except UnicodeDecodeError:
                with open(dir_name, 'rb') as f:
                    results = pickle.load(f, encoding='latin1')
            plots.append(alias)
            self.update_plot('ew', results, color_list[len(plots) - 1], plots)
            if len(plots) ==1:
                self.plot_field(field)
            elif field != self.field_value:
                self.plot_field(field)
                self.field_value = field
            elif self.cleared == 'yes':
                self.plot_field(field)
                self.cleared = 'no'
        else:
            print('File not found, make sure QbC_mgii_v6.py was run with this configuration')
 
 
class MainWindow(QMainWindow):
    """MainWindow class"""     
     
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle('Plots Viewer')#set window title
        self.create_tab_window()
        self.show() #makes instance visible
        self.raise_() #raise instance to top of window stack
        QApplication.setStyle(QStyleFactory.create('Cleanlooks'))
 
    def create_tab_window(self):
        #create tab window
        self.tabs = QTabWidget()
 
        #create tabs
        self.tab1 = PlotWindow()
        self.tab2 = PlotWindow()
        self.tab3 = PlotWindow()
        self.tab4 = PlotWindow()
        self.tab5 = PlotWindowZ()
         
        #add tabs
        self.tabs.addTab(self.tab1, 'Plot 1')
        self.tabs.addTab(self.tab2, 'Plot 2')
        self.tabs.addTab(self.tab3, 'Plot 3')
        self.tabs.addTab(self.tab4, 'Plot 4')
        self.tabs.addTab(self.tab5, 'Plot z')
 
        self.setCentralWidget(self.tabs)
 
        #resize main window
        self.resize(1200, 1100)         
 
def main():
 
    viewer = QApplication(sys.argv)#create new application
    viewer_window = MainWindow()
    # viewer_window.show() #makes instnce visible
    # viewer_window.raise_() #raise instance to top of window stack
    viewer.exec_()#monitor application for events
 
if __name__ == '__main__':
            main()
