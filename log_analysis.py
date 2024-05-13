# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
February 27th, 2024
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931

    ****************************************************************
    * Requirements:                                                *
    *   python 3.7+                                                *
    *                                                              *
    * Dependencies:                                                *
    *   python matplotlib module:                                  *
    *    - pip3 install matplotlib (if pip manager is installed)   *
    *                                                              *
    * Run methods:                                                 *
    *   - GUI (manipulate variables and run from. Default          *
    *          settings set from this script)                      *
    *                                                              *
    * Notes for Anaconda Spyder IDE users:                         *
    *   - If running this from the Anaconda Spyder IDE, before     *
    *     running you will have to turn on the interactive plot    *
    *     if you want the interactive plot to pop-up. You may do   *
    *     this by executing the following command in the console:  *
    *        %matplotlib qt                                        *
    *     Then to turn back on the inline plotting, where the plots*
    *     will appear in the Plots tab execute the following       *
    *     command in the console:                                  *
    *        %matplotlib inline                                    *
    ****************************************************************
"""
##################################################################################################################
# Python int value to range from 0 to 200 to control the GUI zoom. Where it option was added to account for the  #
# different ways individuals have their OS display settings setup. A GUI_zoom = 100, means that default settings #
# are used, whereas a GUI_zoom = 120 means increase GUI size by 20% and GUI_zoom = 80 means decrease GUI by 20%. #
# Examples:                                                                                                      #
#   GUI_zoom = 100 # use default GUI size                                                                        #
#   GUI_zoom = 75  # decrease GUI size by 25%                                                                    #
#                                                                                                                #
# Update use_GUI and GUI_zoom as desired.                                                                        #
##################################################################################################################
GUI_zoom = 100


##################################################################################################################
# Start defining some analysis modes to load basic settings to help speed up certain tedious analysis that are   #
# repetative. Most of the analysis modes are built with the idea of anaylzing LAMMPS thermo log outputs for      #
# process modeling.                                                                                              #
##################################################################################################################
#----------------------------------------------------------------------------------------------------------------#
# blank_mode built by Josh Kemppainen 2/23/2024 to allow the GUI to load with very little predefined settings.   #
# This blank mode provides a basic template of the minimum number of settings that need to be supplied to the    #
# GUI during intialization.                                                                                      #
#----------------------------------------------------------------------------------------------------------------#
blank_mode = {'logfile': 'UPDATE-ME',       # Default logfile name to load during GUI initialization
              'keywords': ['Step', 'Temp'], # Keywords to find logged data (i.e. the columns of data)
              'sections': 'all',            # Sections from logfile to load (all or 1 or 1,2,3 or 1-3 or 1,2,4-6 or ...)
              'xdata': '',                  # Column name in log file for xdata to plot
              'ydata': '',                  # Column name in log file for ydata to plot
              'xlabel':'',                  # xlabel to set on plot (if empty the 'xdata' column name will be used)
              'ylabel':'',                  # ylabel to set on plot (if empty the 'ydata' column name will be used)
              'xscale'  : '',               # multipler value to scale 'xdata' by (float or fraction or empty)
              'yscale'  : '',               # multipler value to scale 'ydata' by (float or fraction or empty)
              'analysis': []}               # list of different analysis to load (if empty nothing will be loaded) 
                                            #   If providing a sublist for an auto-loaded the indexes are as follows:
                                            #     [<method>, <xlo>, <xhi>, <misc>, <name>], where:
                                            #
                                            #      - <method>  is a string and can be 'average' or 'cursor' or
                                            #                  'linear regression' or 'moving average' or 'skip'.
                                            #
                                            #      - <xlo>     is a float to set the lower boundary of the
                                            #                  region to apply the <analysis-type>. 
                                            #
                                            #      - <xhi>     is a float to set the upper boundary of the
                                            #                  region to apply the <analysis-type>. 
                                            #
                                            #       - <misc>   is a string to over ride certain default values,
                                            #                  such as the window of data points to use in the
                                            #                  'moving average'. For more details please run
                                            #                  code and click the "analysis help" button.
                                            #
                                            #       - <name>   is a string to set the name of the analysis in the
                                            #                  legend of the plot.
                                            #   Look at the other modes for examples on how to build a custom
                                            #   analysi mode or to adjust the analysis mode for your models.

#----------------------------------------------------------------------------------------------------------------#
# density_mode built by Josh Kemppainen 2/23/2024 as an example density analysis mode that is loadable. These    #
# settings maybe adjusted as needed. For example the 'xdata' is set as Step and then an 'xscale' value of        #
# '1/2000' is applied to convert the timesteps to picoseconds since a 0.5 fs timestep was used to produce the    #
# example log file.                                                                                              #
#----------------------------------------------------------------------------------------------------------------#
density_mode = {'logfile' : 'EXAMPLES/log_analysis/property=density_ts=0.5.log.lammps',
                'keywords': ['Step', 'Temp'],
                'sections': '1,2',
                'xdata': 'Step',
                'ydata': 'Density',
                'xlabel':'Time (ps)',
                'ylabel':'Density ($g/cm^3$)',
                'xscale'  : '1/2000',
                'yscale'  : '', 
                'analysis': [['average', 1500, 2100, '', 'Density average']]
                }

#----------------------------------------------------------------------------------------------------------------#
# potential_energy_mode built by Josh Kemppainen 2/23/2024 as an example potential energy analysis mode that is  #
# loadable. These settings maybe adjusted as needed. For example the 'xdata' is set as Step and then an 'xscale' #
# value of '1/2000' is applied to convert the timesteps to picoseconds since a 0.5 fs timestep was used to       #
# produce the example log file.                                                                                  #
#----------------------------------------------------------------------------------------------------------------#
potential_energy_mode = {'logfile' : 'EXAMPLES/log_analysis/property=density_ts=0.5.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': '1,2',
                         'xdata': 'Step',
                         'ydata': 'PotEng',
                         'xlabel':'Time (ps)',
                         'ylabel':'Potential Energy (Kcal/mol)',
                         'xscale'  : '1/2000',
                         'yscale'  : '', 
                         'analysis': []
                        }

#----------------------------------------------------------------------------------------------------------------#
# tensile_x_mode built by Josh Kemppainen 2/27/2024 as an example tensile modulus analysis mode that is loadable.#
# In the example below the 'xdata' used for strain data is in the 'v_truex' column and the 'ydata' used for      #
# stress data is in the 'f_sxx_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
tensile_modulus_x_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
                          'keywords': ['Step', 'Temp'],
                          'sections': '1',
                          'xdata': 'v_etruex',
                          'ydata': 'f_sxx_ave',
                          'xlabel':'True Strain',
                          'ylabel':'True Stress (MPa)',
                          'xscale'  : '',
                          'yscale'  : '', 
                          'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                       ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                       ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                         }
#----------------------------------------------------------------------------------------------------------------#
# tensile_x_toughness_mode built by Josh Kemppainen 4/18/2024 as an example tensile analysis mode to compute the #
# toughness (area under the cruve) that is loadable.                                                             #
#----------------------------------------------------------------------------------------------------------------#
tensile_x_toughness_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
                            'keywords': ['Step', 'Temp'],
                            'sections': '1',
                            'xdata': 'v_etruex',
                            'ydata': 'f_sxx_ave',
                            'xlabel':'True Strain',
                            'ylabel':'True Stress (MPa)',
                            'xscale'  : '',
                            'yscale'  : '', 
                            'analysis': [['spline-integration', 0, 0.1, 'window=100; shift=True;', 'Toughness']]
                         }

#----------------------------------------------------------------------------------------------------------------#
# tensile_x_poisson_ratio_ij_mode built by Josh Kemppainen 2/27/2024 as an example poisson's ratio mode that is  #
# loadable. In the example below the 'xdata' used for strain data is in the 'v_truex' column and the 'ydata' is  #
# dependant on if we are going to compute nu_xy or nu_xz ('v_etruey' and 'v_etruez' respectively). The units of  #
# strain is unitless. You may need to update these depending on how you logged your results. Note the usage of   #
# <misc> is empty, thus we are defaulting to not shifting any of the data.                                       #
#----------------------------------------------------------------------------------------------------------------#
tensile_x_poisson_ratio_xy_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruex',
                                   'ydata': 'v_etruey',
                                   'xlabel':'True Strain (in X)',
                                   'ylabel':'True Strain (in Y)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{xy}$  ']]
                                  }
tensile_x_poisson_ratio_xz_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_x_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruex',
                                   'ydata': 'v_etruez',
                                   'xlabel':'True Strain (in X)',
                                   'ylabel':'True Strain (in Z)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{xz}$  ']]
                                  }

#----------------------------------------------------------------------------------------------------------------#
# tensile_y_mode built by Josh Kemppainen 2/27/2024 as an example tensile modulus analysis mode that is loadable.#
# In the example below the 'xdata' used for strain data is in the 'v_truey' column and the 'ydata' used for      #
# stress data is in the 'f_syy_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
tensile_modulus_y_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
                          'keywords': ['Step', 'Temp'],
                          'sections': '1',
                          'xdata': 'v_etruey',
                          'ydata': 'f_syy_ave',
                          'xlabel':'True Strain',
                          'ylabel':'True Stress (MPa)',
                          'xscale'  : '',
                          'yscale'  : '', 
                          'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                       ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                      ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                         }
#----------------------------------------------------------------------------------------------------------------#
# tensile_x_toughness_mode built by Josh Kemppainen 4/18/2024 as an example tensile analysis mode to compute the #
# toughness (area under the cruve) that is loadable.                                                             #
#----------------------------------------------------------------------------------------------------------------#
tensile_y_toughness_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
                            'keywords': ['Step', 'Temp'],
                            'sections': '1',
                            'xdata': 'v_etruey',
                            'ydata': 'f_syy_ave',
                            'xlabel':'True Strain',
                            'ylabel':'True Stress (MPa)',
                            'xscale'  : '',
                            'yscale'  : '', 
                            'analysis': [['spline-integration', 0, 0.1, 'window=100; shift=True;', 'Toughness']]
                         }

#----------------------------------------------------------------------------------------------------------------#
# tensile_y_poisson_ratio_ij_mode built by Josh Kemppainen 2/27/2024 as an example poisson's ratio mode that is  #
# loadable. In the example below the 'xdata' used for strain data is in the 'v_truey' column and the 'ydata' is  #
# dependant on if we are going to compute nu_yx or nu_yz ('v_etruex' and 'v_etruez' respectively). The units of  #
# strain is unitless. You may need to update these depending on how you logged your results. Note the usage of   #
# <misc> is empty, thus we are defaulting to not shifting any of the data.                                       #
#----------------------------------------------------------------------------------------------------------------#
tensile_y_poisson_ratio_yx_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruey',
                                   'ydata': 'v_etruex',
                                   'xlabel':'True Strain (in Y)',
                                   'ylabel':'True Strain (in X)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{yx}$  ']]
                                  }
tensile_y_poisson_ratio_yz_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_y_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruey',
                                   'ydata': 'v_etruez',
                                   'xlabel':'True Strain (in Y)',
                                   'ylabel':'True Strain (in Z)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{yz}$  ']]
                                  }

#----------------------------------------------------------------------------------------------------------------#
# tensile_z_mode built by Josh Kemppainen 2/27/2024 as an example tensile modulus analysis mode that is loadable.#
# In the example below the 'xdata' used for strain data is in the 'v_truez' column and the 'ydata' used for      #
# stress data is in the 'f_szz_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
tensile_modulus_z_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_z_strain_rate=2e8.log.lammps',
                          'keywords': ['Step', 'Temp'],
                          'sections': '1',
                          'xdata': 'v_etruez',
                          'ydata': 'f_szz_ave',
                          'xlabel':'True Strain',
                          'ylabel':'True Stress (MPa)',
                          'xscale'  : '',
                          'yscale'  : '', 
                          'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                       ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                       ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                         }
#----------------------------------------------------------------------------------------------------------------#
# tensile_z_toughness_mode built by Josh Kemppainen 4/18/2024 as an example tensile analysis mode to compute the #
# toughness (area under the cruve) that is loadable.                                                             #
#----------------------------------------------------------------------------------------------------------------#
tensile_z_toughness_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_z_strain_rate=2e8.log.lammps',
                            'keywords': ['Step', 'Temp'],
                            'sections': '1',
                            'xdata': 'v_etruez',
                            'ydata': 'f_szz_ave',
                            'xlabel':'True Strain',
                            'ylabel':'True Stress (MPa)',
                            'xscale'  : '',
                            'yscale'  : '', 
                            'analysis': [['spline-integration', 0, 0.1, 'window=100; shift=True;', 'Toughness']]
                         }

#----------------------------------------------------------------------------------------------------------------#
# tensile_y_poisson_ratio_ij_mode built by Josh Kemppainen 2/27/2024 as an example poisson's ratio mode that is  #
# loadable. In the example below the 'xdata' used for strain data is in the 'v_truez' column and the 'ydata' is  #
# dependant on if we are going to compute nu_zx or nu_zy ('v_etruex' and 'v_etruey' respectively). The units of  #
# strain is unitless. You may need to update these depending on how you logged your results. Note the usage of   #
# <misc> is empty, thus we are defaulting to not shifting any of the data.                                       #
#----------------------------------------------------------------------------------------------------------------#
tensile_z_poisson_ratio_zx_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_z_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruez',
                                   'ydata': 'v_etruex',
                                   'xlabel':'True Strain (in Z)',
                                   'ylabel':'True Strain (in X)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{zx}$  ']]
                                  }
tensile_z_poisson_ratio_zy_mode = {'logfile' : 'EXAMPLES/log_analysis/property=tensile_modulus_z_strain_rate=2e8.log.lammps',
                                   'keywords': ['Step', 'Temp'],
                                   'sections': '1',
                                   'xdata': 'v_etruez',
                                   'ydata': 'v_etruey',
                                   'xlabel':'True Strain (in Z)',
                                   'ylabel':'True Strain (in Y)',
                                   'xscale'  : '',
                                   'yscale'  : '', 
                                   'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                                ['linear regression', 0, 0.02, '', '$\\nu_{zy}$  ']]
                                  }

#----------------------------------------------------------------------------------------------------------------#
# shear_xy_mode built by Josh Kemppainen 2/27/2024 as an example shear modulus analysis mode that is loadable.   #
# In the example below the 'xdata' used for strain data is in the 'v_etruexy' column and the 'ydata' used for    #
# stress data is in the 'f_sxy_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
shear_modulus_xy_mode = {'logfile' : 'EXAMPLES/log_analysis/property=shear_modulus_xy_strain_rate=2e8.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': '1',
                         'xdata': 'v_etruexy',
                         'ydata': 'f_sxy_ave',
                         'xlabel':'True Strain',
                         'ylabel':'True Stress (MPa)',
                         'xscale'  : '',
                         'yscale'  : '', 
                         'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                      ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                      ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                        }

#----------------------------------------------------------------------------------------------------------------#
# shear_xy_mode built by Josh Kemppainen 2/27/2024 as an example shear modulus analysis mode that is loadable.   #
# In the example below the 'xdata' used for strain data is in the 'v_etruexy' column and the 'ydata' used for    #
# stress data is in the 'f_sxy_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
shear_modulus_xz_mode = {'logfile' : 'EXAMPLES/log_analysis/property=shear_modulus_xz_strain_rate=2e8.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': '1',
                         'xdata': 'v_etruexz',
                         'ydata': 'f_sxz_ave',
                         'xlabel':'True Strain',
                         'ylabel':'True Stress (MPa)',
                         'xscale'  : '',
                         'yscale'  : '', 
                         'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                      ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                      ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                        }

#----------------------------------------------------------------------------------------------------------------#
# shear_yz_mode built by Josh Kemppainen 2/27/2024 as an example shear modulus analysis mode that is loadable.   #
# In the example below the 'xdata' used for strain data is in the 'v_etrueyz' column and the 'ydata' used for    #
# stress data is in the 'f_syz_ave' column. The units of the stress data is in MPa and strain is unitless. You   #
# may need to update these depending on how you logged your results. Note the usage of <misc> = 'shift=True'     #
# which is used to shift all the y-data by the y-intercept of the linear regression model. This is to allow for  #
# the prediction of yield strength, where the removal of residual stresses is being performed.                   #
#----------------------------------------------------------------------------------------------------------------#
shear_modulus_yz_mode = {'logfile' : 'EXAMPLES/log_analysis/property=shear_modulus_yz_strain_rate=2e8.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': '1',
                         'xdata': 'v_etrueyz',
                         'ydata': 'f_syz_ave',
                         'xlabel':'True Strain',
                         'ylabel':'True Stress (MPa)',
                         'xscale'  : '',
                         'yscale'  : '', 
                         'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                      ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                      ['cursor', '', '', 'x=0.02; y=5;', 'Yield strength']]
                        }

#----------------------------------------------------------------------------------------------------------------#
# This is special mode as it will show how to use strings to perform vectorized manipulation of arrays to compute#
# a new array of data to plot. In general the 'xdata' or 'ydata' keys should be columns in the logfile, however  #
# users may define a compute by setting a string as 'compute: MATH', where MATH is math that the Python eval()   #
# function can handle. A little bit of understanding needs to be understood about how things are setup in the    #
# code to make use of this option detailed below:                                                                #
#   - data is setup in a dictionary as d = {'COLUMN-name': numpy.array(data from logfile)}                       #
#   - numpy can perform vectorized math operations where the + or - or * or / or ** operations can be used       #
#                                                                                                                #
# As an example say a LAMMPS logfile had columns x and y, with x-values = [1,2,3] and y-values = [4,5,6]. The d  #
# dictionary is setup as d = {'x': numpy.array([1,2,3]), 'y': numpy.array([4,5,6])}. Where any math operation    #
# that is compatable with the eval() command and numpy vectorization can be used. You can google around for to   #
# find this, but in general almost any basic math opertation is supported. Here are a few examples of setting    #
# compute strings:                                                                                               #
#    "compute: d['x']    + d['y']"        the code would compute a new array = [5,7,9].                          #
#    "compute: d['x']**2 + d['y']"        the code would compute a new array = [5,9,15].                         #
#    "compute: d['x']**2 / d['y']**2"     the code would compute a new array = [0.0625, 0.16, 0.25].             #
#    "compute: d['x']**2 * d['y']**0.5"   the code would compute a new array = [2.0, 8.94, 22.05].               #
#                                                                                                                #
# You may then set the 'xdata' and/or 'ydata' in the mode dictionaries as these "compute" strings, which will    #
# tell the code to perform those math operations, using the data set by the columns of the logfile.              #
#                                                                                                                #
# This can be useful to compute a new array of data, where the most obvious case is the need to compute the von  #
# mises stress criteria from a shear simulation. As in some research groups shear simulations are performed and  #
# the equivalent tensile yield strength is computed from the von mises stress criteria. Below are some notes     #
# on the von mises stress criteria equation and how the example shear stress log file is setup.                  #
#    svm = sigma_von_mises (computed)                                                                            #
#    s11 = sigma_11 (column in logfile: f_sxx_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s22 = sigma_22 (column in logfile: f_syy_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s33 = sigma_33 (column in logfile: f_szz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s12 = sigma_12 (column in logfile: f_sxy_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s23 = sigma_23 (column in logfile: f_syz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
#    s31 = sigma_31 (column in logfile: f_sxz_ave  PLEASE NOTE THIS MAY NOT BE THE SAME IN YOUR LOGFILES)        #
# svm = ( ((s11-s22)**2 + (s22-s33)**2 + (s33-s11)**2 + 6*(s12**2 + s23**2 + s31**2))/2 )**0.5                   #
#----------------------------------------------------------------------------------------------------------------#
# Method 1: Use a "brute force" typing of string (this serves as an example)
svm = "( ((d['f_sxx_ave']-d['f_syy_ave'])**2 + (d['f_syy_ave']-d['f_szz_ave'])**2 + (d['f_szz_ave']-d['f_sxx_ave'])**2 + 6*(d['f_sxy_ave']**2 + d['f_syz_ave']**2 + d['f_sxz_ave']**2))/2 )**0.5"

# Method 2: Using formated strings, this method is prefered as it lets others change column names a bit easier
s11 = "d['{}']".format('f_sxx_ave')
s22 = "d['{}']".format('f_syy_ave')
s33 = "d['{}']".format('f_szz_ave')
s12 = "d['{}']".format('f_sxy_ave')
s23 = "d['{}']".format('f_syz_ave')
s31 = "d['{}']".format('f_sxz_ave')
svm = "( (({s11}-{s22})**2 + ({s22}-{s33})**2 + ({s33}-{s11})**2 + 6*({s23}**2 + {s31}**2 + {s12}**2))/2 )**0.5".format(s11=s11, s22=s22, s33=s33, s12=s12, s23=s23, s31=s31)

# Define the von mises mode
shear_von_mises_mode = {'logfile' : 'EXAMPLES/log_analysis/property=shear_modulus_xz_strain_rate=2e8.log.lammps',
                        'keywords': ['Step', 'Temp'],
                        'sections': '1',
                        'xdata': 'v_etruexz',
                        'ydata': f'compute: {svm}', # use variable from above since this is a long equation
                        'xlabel':'True Strain',
                        'ylabel':'Von Mises Stress (MPa)',
                        'xscale'  : '',
                        'yscale'  : '', 
                        'analysis': [['moving average', 0, 0.1, 'window=50', 'moving average'],
                                     ['linear regression', 0, 0.02, 'shift=True', 'Modulus'],
                                     ['cursor', '', '', 'x=0; y=0;', 'Yield strength']]
                        }

#----------------------------------------------------------------------------------------------------------------#
# bulk_modulus_v0_v1_mode built by Josh Kemppainen 2/27/2024 as an example density analysis mode that is         #
# loadable. These settings maybe adjusted as needed.                                                             #
#----------------------------------------------------------------------------------------------------------------#
bulk_modulus_v0_v1_mode = {'logfile' : 'EXAMPLES\log_analysis\property=bulk_modulus_v0_v1_calc_ts=0.5_p0=1atm_p1=5000atm.log.lammps',
                           'keywords': ['Step', 'Temp'],
                           'sections': '1-3',
                           'xdata': 'Step',
                           'ydata': 'Volume',
                           'xlabel':'Time (ps)',
                           'ylabel':'Volume (A$^3$)',
                           'xscale'  : '1/2000',
                           'yscale'  : '', 
                           'analysis': [['average', 0, 400, '', 'v0 average (A$^3$)'],
                                        ['average', 1000, 2100, '', 'v1 average (A$^3$)']]
                            }


#----------------------------------------------------------------------------------------------------------------#
# Tg_CTE_regression_mode built by Josh Kemppainen 2/27/2024 as an example density analysis mode that is loadable.#
# These settings maybe adjusted as needed.                                                                       #
#----------------------------------------------------------------------------------------------------------------#
Tg_CTE_regression_mode = {'logfile' : 'EXAMPLES\log_analysis\properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps',
                          'keywords': ['Step', 'Temp'],
                          'sections': 'all',
                          'xdata': 'Temp',
                          'ydata': 'Volume',
                          'xlabel':'Temperature (K)',
                          'ylabel':'Volume (A$^3$)',
                          'xscale'  : '',
                          'yscale'  : '', 
                          'analysis': [['moving average', 100, 800, 'window=100', 'moving average'],
                                       ['average', 290, 310, '', 'room temp volume average (A$^3$)'],
                                       ['linear regression', 100, 300, 'extend=250', 'CTE below Tg'], 
                                       ['linear regression', 600, 750, 'extend=-100', 'CTE above Tg'],
                                       ['cursor', '', '', 'x=460; y=190_000;', 'Tg']]
                            }

#----------------------------------------------------------------------------------------------------------------#
# Tg_CTE_hyperbola_mode built by Josh Kemppainen 4/15/2024 as an example density analysis mode that is loadable. #
# These settings maybe adjusted as needed.                                                                       #
#----------------------------------------------------------------------------------------------------------------#
Tg_CTE_hyperbola_mode = {'logfile' : 'EXAMPLES\log_analysis\properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': 'all',
                         'xdata': 'Temp',
                         'ydata': 'Volume',
                         'xlabel':'Temperature (K)',
                         'ylabel':'Volume (A$^3$)',
                         'xscale'  : '',
                         'yscale'  : '', 
                         'analysis': [['average', 290, 310, '', 'room temp volume average (A$^3$)'],
                                      ['moving average', 100, 800, 'window=100', 'moving average'],
                                      ['hyperbola', 100, 800, 'p=0.9', 'Hyperbola fit']]
                            }

#----------------------------------------------------------------------------------------------------------------#
# Tg_CTE_piecewise_mode built by Josh Kemppainen 4/15/2024 as an example density analysis mode that is loadable.  #
# These settings maybe adjusted as needed.                                                                       #
#----------------------------------------------------------------------------------------------------------------#
Tg_CTE_piecewise_mode = {'logfile' : 'EXAMPLES\log_analysis\properties=Tg_and_CTE_heating_rate=50_k_ns.log.lammps',
                         'keywords': ['Step', 'Temp'],
                         'sections': 'all',
                         'xdata': 'Temp',
                         'ydata': 'Volume',
                         'xlabel':'Temperature (K)',
                         'ylabel':'Volume (A$^3$)',
                         'xscale'  : '',
                         'yscale'  : '', 
                         'analysis': [['average', 290, 310, '', 'room temp volume average (A$^3$)'],
                                      ['moving average', 100, 800, 'window=100', 'moving average'],
                                      ['piecewise-regression', 100, 800, 'n=1', 'Segmented piecewise-regression']]
                            }


#----------------------------------------------------------------------------------------------------------------#
# min_density_mode built by Josh Kemppainen 5/13/2024 as an example density analysis mode that is loadable. These#
# settings maybe adjusted as needed. For example the 'xdata' is set as Step and then an 'xscale' value of        #
# '1/2000' is applied to convert the timesteps to picoseconds since a 0.5 fs timestep was used to produce the    #
# example log file.                                                                                              #
#----------------------------------------------------------------------------------------------------------------#
min_density_mode = {'logfile' : 'EXAMPLES/log_analysis/property=density_ts=0.5.log.lammps',
                    'keywords': ['Step', 'Temp'],
                    'sections': '1,2',
                       'xdata': 'Step',
                       'ydata': 'Density',
                      'xlabel':'Time (ps)',
                      'ylabel':'Density ($g/cm^3$)',
                    'xscale'  : '1/2000',
                    'yscale'  : '', 
                    'analysis': [['minimum', 0, 2100, 'window=10', 'Density minimum']]
                   }

#----------------------------------------------------------------------------------------------------------------#
# min_density_mode built by Josh Kemppainen 5/13/2024 as an example density analysis mode that is loadable. These#
# settings maybe adjusted as needed. For example the 'xdata' is set as Step and then an 'xscale' value of        #
# '1/2000' is applied to convert the timesteps to picoseconds since a 0.5 fs timestep was used to produce the    #
# example log file.                                                                                              #
#----------------------------------------------------------------------------------------------------------------#
max_density_mode = {'logfile' : 'EXAMPLES/log_analysis/property=density_ts=0.5.log.lammps',
                    'keywords': ['Step', 'Temp'],
                    'sections': '1,2',
                       'xdata': 'Step',
                       'ydata': 'Density',
                      'xlabel':'Time (ps)',
                      'ylabel':'Density ($g/cm^3$)',
                    'xscale'  : '1/2000',
                    'yscale'  : '', 
                    'analysis': [['maximum', 0, 2100, '', 'Density maximum']]
                   }

#----------------------------------------------------------------------------------------------------------------#
# Group all of the modes into one dictionary. You may then open each mode in the GUI based on the name provided  #
# as the key to the mode dictionary. If more modes are defined above, you MUST set them in this dictionary and   #
# provide a key to access the mode from within the GUI, when the GUI is running.                                 #                                                                            
#----------------------------------------------------------------------------------------------------------------#
modes = {'blank': blank_mode,
         'density': density_mode,
         'minimum density': min_density_mode,
         'maximum density': max_density_mode,
         'potential energy': potential_energy_mode,
         'tensile modulus x': tensile_modulus_x_mode,
         'tensile modulus y': tensile_modulus_y_mode,
         'tensile modulus z': tensile_modulus_z_mode,
         'tensile x tougness': tensile_x_toughness_mode,
         'tensile y tougness': tensile_y_toughness_mode,
         'tensile z tougness': tensile_z_toughness_mode,
         'poisson ratio xy (tensile x)': tensile_x_poisson_ratio_xy_mode,
         'poisson ratio xz (tensile x)': tensile_x_poisson_ratio_xz_mode,    
         'poisson ratio yx (tensile y)': tensile_y_poisson_ratio_yx_mode,
         'poisson ratio yz (tensile y)': tensile_y_poisson_ratio_yz_mode,
         'poisson ratio zx (tensile z)': tensile_z_poisson_ratio_zx_mode,
         'poisson ratio zy (tensile z)': tensile_z_poisson_ratio_zy_mode,
         'shear modulus xy': shear_modulus_xy_mode,
         'shear modulus xz': shear_modulus_xz_mode,
         'shear modulus yz': shear_modulus_yz_mode,
         'shear von mises': shear_von_mises_mode,
         'bulk modulus v0 and v1': bulk_modulus_v0_v1_mode,
         'Tg and CTE linear regression': Tg_CTE_regression_mode,
         'Tg and CTE hyperbola': Tg_CTE_hyperbola_mode,
         'Tg and CTE piecewise':Tg_CTE_piecewise_mode,
        }


##################################################################################################################
# Define basic settings, such as the 'mode' to load with, the predefined modes, and the figure image quality.    #
#                                                                                                                #
# Update settings as desired.                                                                                    #
##################################################################################################################
settings = {'mode': 'blank',    # Default run mode to load GUI with (could be updated to any of the modes from above)
            'modes': modes,     # Define all possible modes that could be loaded in the GUI
            'save-fig' : True,  # True or False option to save figure are .jpeg
            'image-dpi': 300,   # Set the image quality of the generated image (NOTE higher dpi's will result in slower code)
           }


###################################
### Import needed files and run ###
###################################
if __name__ == "__main__":  
    import src.log_analysis.GUI as GUI
    GUI.GUI(settings, GUI_zoom)