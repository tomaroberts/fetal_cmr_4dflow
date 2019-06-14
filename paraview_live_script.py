#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
cine_vol_masked_t = LegacyVTKReader(FileNames=['C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t0.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t1.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t2.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t3.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t4.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t5.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t6.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t7.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t8.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t9.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t10.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t11.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t12.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t13.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t14.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t15.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t16.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t17.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t18.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t19.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t20.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t21.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t22.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t23.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\cine_vol_masked_t24.vtk'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Legacy VTK Reader'
vel_vol_masked_VxVyVz_t = LegacyVTKReader(FileNames=['C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t0.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t1.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t2.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t3.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t4.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t5.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t6.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t7.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t8.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t9.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t10.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t11.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t12.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t13.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t14.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t15.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t16.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t17.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t18.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t19.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t20.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t21.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t22.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t23.vtk', 'C:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\fcmr194\\vel_vol\\paraview\\vel_vol_masked_VxVy-Vz_t24.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [955, 778]

# show data in view
cine_vol_masked_tDisplay = Show(cine_vol_masked_t, renderView1)
# trace defaults for the display properties.
cine_vol_masked_tDisplay.Representation = 'Outline'
cine_vol_masked_tDisplay.AmbientColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.ColorArrayName = ['POINTS', '']
cine_vol_masked_tDisplay.OSPRayScaleArray = 'magnitude_intensity'
cine_vol_masked_tDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
cine_vol_masked_tDisplay.SelectOrientationVectors = 'None'
cine_vol_masked_tDisplay.ScaleFactor = 7.2
cine_vol_masked_tDisplay.SelectScaleArray = 'magnitude_intensity'
cine_vol_masked_tDisplay.GlyphType = 'Arrow'
cine_vol_masked_tDisplay.GlyphTableIndexArray = 'magnitude_intensity'
cine_vol_masked_tDisplay.DataAxesGrid = 'GridAxesRepresentation'
cine_vol_masked_tDisplay.PolarAxes = 'PolarAxesRepresentation'
cine_vol_masked_tDisplay.ScalarOpacityUnitDistance = 1.7383632950303343
#cine_vol_masked_tDisplay.InputVectors = [None, '']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
cine_vol_masked_tDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
cine_vol_masked_tDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
cine_vol_masked_tDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# show data in view
vel_vol_masked_VxVyVz_tDisplay = Show(vel_vol_masked_VxVyVz_t, renderView1)
# trace defaults for the display properties.
vel_vol_masked_VxVyVz_tDisplay.Representation = 'Outline'
vel_vol_masked_VxVyVz_tDisplay.AmbientColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.ColorArrayName = ['POINTS', '']
vel_vol_masked_VxVyVz_tDisplay.OSPRayScaleArray = 'velocity_magnitude'
vel_vol_masked_VxVyVz_tDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
vel_vol_masked_VxVyVz_tDisplay.SelectOrientationVectors = 'vector_field'
vel_vol_masked_VxVyVz_tDisplay.ScaleFactor = 7.2
vel_vol_masked_VxVyVz_tDisplay.SelectScaleArray = 'velocity_magnitude'
vel_vol_masked_VxVyVz_tDisplay.GlyphType = 'Arrow'
vel_vol_masked_VxVyVz_tDisplay.GlyphTableIndexArray = 'velocity_magnitude'
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid = 'GridAxesRepresentation'
vel_vol_masked_VxVyVz_tDisplay.PolarAxes = 'PolarAxesRepresentation'
vel_vol_masked_VxVyVz_tDisplay.ScalarOpacityUnitDistance = 1.7383632950303343
#vel_vol_masked_VxVyVz_tDisplay.InputVectors = ['POINTS', 'vector_field']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
vel_vol_masked_VxVyVz_tDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
vel_vol_masked_VxVyVz_tDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(cine_vol_masked_t)

# set scalar coloring
ColorBy(cine_vol_masked_tDisplay, ('POINTS', 'magnitude_intensity'))

# rescale color and/or opacity maps used to include current data range
cine_vol_masked_tDisplay.RescaleTransferFunctionToDataRange(True, True)

# change representation type
cine_vol_masked_tDisplay.SetRepresentationType('Volume')

# get color transfer function/color map for 'magnitude_intensity'
magnitude_intensityLUT = GetColorTransferFunction('magnitude_intensity')
magnitude_intensityLUT.RGBPoints = [-7.28000020980835, 0.231373, 0.298039, 0.752941, 87.65499806404114, 0.865003, 0.865003, 0.865003, 182.58999633789062, 0.705882, 0.0156863, 0.14902]
magnitude_intensityLUT.ScalarRangeInitialized = 1.0

# Properties modified on cine_vol_masked_tDisplay
cine_vol_masked_tDisplay.SelectMapper = 'Resample To Image'

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
magnitude_intensityLUT.ApplyPreset('Grayscale', True)

# set active source
SetActiveSource(vel_vol_masked_VxVyVz_t)

# create a new 'Glyph'
glyph1 = Glyph(Input=vel_vol_masked_VxVyVz_t,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'velocity_magnitude']
glyph1.Vectors = ['POINTS', 'vector_field']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 7.2
glyph1.GlyphMode = 'All Points'
glyph1.MaximumNumberOfSamplePoints = 10000000
glyph1.GlyphTransform = 'Transform2'

# hide data in view
Hide(vel_vol_masked_VxVyVz_t, renderView1)

# Properties modified on glyph1
glyph1.ScaleFactor = 0.1

# get color transfer function/color map for 'velocity_magnitude'
velocity_magnitudeLUT = GetColorTransferFunction('velocity_magnitude')
velocity_magnitudeLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 11.645000457763672, 0.865003, 0.865003, 0.865003, 23.290000915527344, 0.705882, 0.0156863, 0.14902]
velocity_magnitudeLUT.ScalarRangeInitialized = 1.0

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'velocity_magnitude']
glyph1Display.LookupTable = velocity_magnitudeLUT
glyph1Display.OSPRayScaleArray = 'velocity_magnitude'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 7.2
glyph1Display.SelectScaleArray = 'velocity_magnitude'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'velocity_magnitude'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 3.6
glyph1Display.SetScaleArray = ['POINTS', 'velocity_magnitude']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'velocity_magnitude']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
#glyph1Display.InputVectors = ['POINTS', 'GlyphVector']

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# Rescale transfer function
velocity_magnitudeLUT.RescaleTransferFunction(0.0, 30.0)

# get opacity transfer function/opacity map for 'velocity_magnitude'
velocity_magnitudePWF = GetOpacityTransferFunction('velocity_magnitude')
velocity_magnitudePWF.Points = [0.0, 0.0, 0.5, 0.0, 23.290000915527344, 1.0, 0.5, 0.0]
velocity_magnitudePWF.ScalarRangeInitialized = 1

# Rescale transfer function
velocity_magnitudePWF.RescaleTransferFunction(0.0, 30.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
velocity_magnitudeLUT.ApplyPreset('Rainbow Desaturated', True)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [955, 778]

# current camera placement for renderView1
renderView1.CameraPosition = [47.25540552915633, 105.23644996989677, 28.34453809869176]
renderView1.CameraFocalPoint = [-3.126612682896695, -111.29896476971966, 32.81751852599477]
renderView1.CameraViewUp = [0.16504597062440005, -0.058739378433129194, -0.9845351761120267]
renderView1.CameraParallelScale = 57.5521502639128

# save screenshot
SaveScreenshot('C:/Users/tr17/Documents/Projects/PC_Fetal_CMR/Data/fcmr194/vel_vol/paraview/test.png', renderView1, ImageResolution=[955, 778])

# save state
SaveState('C:/Users/tr17/Documents/Projects/PC_Fetal_CMR/Data/fcmr194/vel_vol/paraview/test_state.pvsm')