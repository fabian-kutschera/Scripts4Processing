# Project 3D Surface onto 2D plane following
# https://discourse.paraview.org/t/project-3d-surface-onto-2d-plane/709
# https://stackoverflow.com/questions/45544709/paraview-view-a-2d-projection-of-a-3d-object/45595967#45595967
# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
hFFtestsurfacexdmf = XDMFReader(registrationName='HFFtest-surface.xdmf', FileNames=['/import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/simple_fault_geometry/Simple_East_M7.341/HFFtest-surface.xdmf'])
hFFtestsurfacexdmf.CellArrayStatus = ['U', 'V', 'W', 'partition', 'u', 'v', 'w']
hFFtestsurfacexdmf.GridStatus = ['step_000000000000', 'step_000000000001', 'step_000000000002', 'step_000000000003', 'step_000000000004', 'step_000000000005', 'step_000000000006', 'step_000000000007', 'step_000000000008', 'step_000000000009', 'step_000000000010', 'step_000000000011', 'step_000000000012', 'step_000000000013', 'step_000000000014', 'step_000000000015', 'step_000000000016', 'step_000000000017', 'step_000000000018', 'step_000000000019', 'step_000000000020', 'step_000000000021', 'step_000000000022', 'step_000000000023', 'step_000000000024', 'step_000000000025', 'step_000000000026', 'step_000000000027', 'step_000000000028', 'step_000000000029', 'step_000000000030', 'step_000000000031', 'step_000000000032', 'step_000000000033', 'step_000000000034', 'step_000000000035', 'step_000000000036', 'step_000000000037', 'step_000000000038', 'step_000000000039', 'step_000000000040', 'step_000000000041', 'step_000000000042', 'step_000000000043', 'step_000000000044', 'step_000000000045', 'step_000000000046', 'step_000000000047', 'step_000000000048', 'step_000000000049', 'step_000000000050', 'step_000000000051', 'step_000000000052', 'step_000000000053', 'step_000000000054', 'step_000000000055', 'step_000000000056', 'step_000000000057', 'step_000000000058', 'step_000000000059', 'step_000000000060', 'step_000000000061', 'step_000000000062', 'step_000000000063', 'step_000000000064', 'step_000000000065', 'step_000000000066', 'step_000000000067', 'step_000000000068', 'step_000000000069', 'step_000000000070', 'step_000000000071', 'step_000000000072', 'step_000000000073', 'step_000000000074', 'step_000000000075', 'step_000000000076', 'step_000000000077', 'step_000000000078', 'step_000000000079', 'step_000000000080', 'step_000000000081', 'step_000000000082', 'step_000000000083', 'step_000000000084', 'step_000000000085', 'step_000000000086', 'step_000000000087', 'step_000000000088', 'step_000000000089', 'step_000000000090', 'step_000000000091', 'step_000000000092', 'step_000000000093', 'step_000000000094', 'step_000000000095', 'step_000000000096', 'step_000000000097', 'step_000000000098', 'step_000000000099', 'step_000000000100', 'step_000000000101', 'step_000000000102', 'step_000000000103', 'step_000000000104', 'step_000000000105', 'step_000000000106', 'step_000000000107', 'step_000000000108', 'step_000000000109', 'step_000000000110', 'step_000000000111', 'step_000000000112', 'step_000000000113', 'step_000000000114', 'step_000000000115', 'step_000000000116', 'step_000000000117', 'step_000000000118', 'step_000000000119', 'step_000000000120', 'step_000000000121', 'step_000000000122', 'step_000000000123', 'step_000000000124', 'step_000000000125', 'step_000000000126', 'step_000000000127', 'step_000000000128', 'step_000000000129', 'step_000000000130', 'step_000000000131', 'step_000000000132', 'step_000000000133', 'step_000000000134', 'step_000000000135', 'step_000000000136', 'step_000000000137', 'step_000000000138', 'step_000000000139', 'step_000000000140', 'step_000000000141', 'step_000000000142', 'step_000000000143', 'step_000000000144', 'step_000000000145', 'step_000000000146', 'step_000000000147', 'step_000000000148', 'step_000000000149', 'step_000000000150', 'step_000000000151', 'step_000000000152', 'step_000000000153', 'step_000000000154', 'step_000000000155', 'step_000000000156', 'step_000000000157', 'step_000000000158', 'step_000000000159', 'step_000000000160', 'step_000000000161', 'step_000000000162', 'step_000000000163', 'step_000000000164', 'step_000000000165', 'step_000000000166', 'step_000000000167', 'step_000000000168', 'step_000000000169', 'step_000000000170', 'step_000000000171', 'step_000000000172', 'step_000000000173', 'step_000000000174', 'step_000000000175', 'step_000000000176', 'step_000000000177', 'step_000000000178', 'step_000000000179', 'step_000000000180', 'step_000000000181', 'step_000000000182', 'step_000000000183', 'step_000000000184', 'step_000000000185', 'step_000000000186', 'step_000000000187', 'step_000000000188', 'step_000000000189', 'step_000000000190', 'step_000000000191', 'step_000000000192', 'step_000000000193', 'step_000000000194', 'step_000000000195', 'step_000000000196', 'step_000000000197', 'step_000000000198', 'step_000000000199', 'step_000000000200', 'step_000000000201', 'step_000000000202', 'step_000000000203', 'step_000000000204', 'step_000000000205', 'step_000000000206', 'step_000000000207', 'step_000000000208', 'step_000000000209', 'step_000000000210', 'step_000000000211', 'step_000000000212', 'step_000000000213', 'step_000000000214', 'step_000000000215', 'step_000000000216', 'step_000000000217', 'step_000000000218', 'step_000000000219', 'step_000000000220', 'step_000000000221', 'step_000000000222', 'step_000000000223', 'step_000000000224', 'step_000000000225', 'step_000000000226', 'step_000000000227', 'step_000000000228', 'step_000000000229', 'step_000000000230', 'step_000000000231', 'step_000000000232', 'step_000000000233', 'step_000000000234', 'step_000000000235', 'step_000000000236', 'step_000000000237', 'step_000000000238', 'step_000000000239', 'step_000000000240', 'step_000000000241', 'step_000000000242', 'step_000000000243', 'step_000000000244', 'step_000000000245', 'step_000000000246', 'step_000000000247', 'step_000000000248', 'step_000000000249', 'step_000000000250', 'step_000000000251', 'step_000000000252', 'step_000000000253', 'step_000000000254', 'step_000000000255', 'step_000000000256', 'step_000000000257', 'step_000000000258', 'step_000000000259', 'step_000000000260', 'step_000000000261', 'step_000000000262', 'step_000000000263', 'step_000000000264', 'step_000000000265', 'step_000000000266', 'step_000000000267', 'step_000000000268', 'step_000000000269', 'step_000000000270', 'step_000000000271', 'step_000000000272', 'step_000000000273', 'step_000000000274', 'step_000000000275', 'step_000000000276', 'step_000000000277', 'step_000000000278', 'step_000000000279', 'step_000000000280', 'step_000000000281', 'step_000000000282', 'step_000000000283', 'step_000000000284', 'step_000000000285', 'step_000000000286', 'step_000000000287', 'step_000000000288', 'step_000000000289', 'step_000000000290', 'step_000000000291', 'step_000000000292', 'step_000000000293', 'step_000000000294', 'step_000000000295', 'step_000000000296', 'step_000000000297', 'step_000000000298', 'step_000000000299', 'step_000000000300', 'step_000000000301', 'step_000000000302', 'step_000000000303', 'step_000000000304', 'step_000000000305', 'step_000000000306', 'step_000000000307', 'step_000000000308', 'step_000000000309', 'step_000000000310', 'step_000000000311', 'step_000000000312', 'step_000000000313', 'step_000000000314', 'step_000000000315', 'step_000000000316', 'step_000000000317', 'step_000000000318', 'step_000000000319', 'step_000000000320', 'step_000000000321', 'step_000000000322', 'step_000000000323', 'step_000000000324', 'step_000000000325', 'step_000000000326', 'step_000000000327', 'step_000000000328', 'step_000000000329', 'step_000000000330', 'step_000000000331', 'step_000000000332', 'step_000000000333', 'step_000000000334', 'step_000000000335', 'step_000000000336', 'step_000000000337', 'step_000000000338', 'step_000000000339', 'step_000000000340', 'step_000000000341', 'step_000000000342', 'step_000000000343', 'step_000000000344', 'step_000000000345', 'step_000000000346', 'step_000000000347', 'step_000000000348', 'step_000000000349', 'step_000000000350', 'step_000000000351', 'step_000000000352', 'step_000000000353', 'step_000000000354', 'step_000000000355', 'step_000000000356', 'step_000000000357', 'step_000000000358', 'step_000000000359', 'step_000000000360', 'step_000000000361', 'step_000000000362', 'step_000000000363', 'step_000000000364', 'step_000000000365', 'step_000000000366', 'step_000000000367', 'step_000000000368', 'step_000000000369', 'step_000000000370', 'step_000000000371', 'step_000000000372', 'step_000000000373', 'step_000000000374', 'step_000000000375', 'step_000000000376', 'step_000000000377', 'step_000000000378', 'step_000000000379', 'step_000000000380', 'step_000000000381', 'step_000000000382', 'step_000000000383', 'step_000000000384', 'step_000000000385', 'step_000000000386', 'step_000000000387', 'step_000000000388', 'step_000000000389', 'step_000000000390', 'step_000000000391', 'step_000000000392', 'step_000000000393', 'step_000000000394', 'step_000000000395', 'step_000000000396', 'step_000000000397', 'step_000000000398', 'step_000000000399', 'step_000000000400']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
hFFtestsurfacexdmfDisplay = Show(hFFtestsurfacexdmf, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'partition'
partitionLUT = GetColorTransferFunction('partition')

# get opacity transfer function/opacity map for 'partition'
partitionPWF = GetOpacityTransferFunction('partition')

# trace defaults for the display properties.
hFFtestsurfacexdmfDisplay.Representation = 'Surface'
hFFtestsurfacexdmfDisplay.ColorArrayName = ['CELLS', 'partition']
hFFtestsurfacexdmfDisplay.LookupTable = partitionLUT
hFFtestsurfacexdmfDisplay.SelectTCoordArray = 'None'
hFFtestsurfacexdmfDisplay.SelectNormalArray = 'None'
hFFtestsurfacexdmfDisplay.SelectTangentArray = 'None'
hFFtestsurfacexdmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
hFFtestsurfacexdmfDisplay.SelectOrientationVectors = 'None'
hFFtestsurfacexdmfDisplay.ScaleFactor = 20500.0
hFFtestsurfacexdmfDisplay.SelectScaleArray = 'partition'
hFFtestsurfacexdmfDisplay.GlyphType = 'Arrow'
hFFtestsurfacexdmfDisplay.GlyphTableIndexArray = 'partition'
hFFtestsurfacexdmfDisplay.GaussianRadius = 1025.0
hFFtestsurfacexdmfDisplay.SetScaleArray = [None, '']
hFFtestsurfacexdmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
hFFtestsurfacexdmfDisplay.OpacityArray = [None, '']
hFFtestsurfacexdmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
hFFtestsurfacexdmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
hFFtestsurfacexdmfDisplay.PolarAxes = 'PolarAxesRepresentation'
hFFtestsurfacexdmfDisplay.ScalarOpacityFunction = partitionPWF
hFFtestsurfacexdmfDisplay.ScalarOpacityUnitDistance = 3530.290776620291
hFFtestsurfacexdmfDisplay.OpacityArrayName = ['CELLS', 'partition']

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
hFFtestsurfacexdmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToLast()

# create a new 'NetCDF Reader'
tanioka_simple_east_M7343nc = NetCDFReader(registrationName='tanioka_simple_east_M7.343.nc', FileName=['/import/freenas-m-05-seissol/kutschera/HIWI/samoa/tanioka_simple_east_M7.343.nc'])
tanioka_simple_east_M7343nc.Dimensions = '(y, x)'

# show data in view
tanioka_simple_east_M7343ncDisplay = Show(tanioka_simple_east_M7343nc, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'z'
zLUT = GetColorTransferFunction('z')

# get opacity transfer function/opacity map for 'z'
zPWF = GetOpacityTransferFunction('z')

# trace defaults for the display properties.
tanioka_simple_east_M7343ncDisplay.Representation = 'Slice'
tanioka_simple_east_M7343ncDisplay.ColorArrayName = ['POINTS', 'z']
tanioka_simple_east_M7343ncDisplay.LookupTable = zLUT
tanioka_simple_east_M7343ncDisplay.SelectTCoordArray = 'None'
tanioka_simple_east_M7343ncDisplay.SelectNormalArray = 'None'
tanioka_simple_east_M7343ncDisplay.SelectTangentArray = 'None'
tanioka_simple_east_M7343ncDisplay.OSPRayScaleArray = 'z'
tanioka_simple_east_M7343ncDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
tanioka_simple_east_M7343ncDisplay.SelectOrientationVectors = 'None'
tanioka_simple_east_M7343ncDisplay.ScaleFactor = 20500.0
tanioka_simple_east_M7343ncDisplay.SelectScaleArray = 'None'
tanioka_simple_east_M7343ncDisplay.GlyphType = 'Arrow'
tanioka_simple_east_M7343ncDisplay.GlyphTableIndexArray = 'None'
tanioka_simple_east_M7343ncDisplay.GaussianRadius = 1025.0
tanioka_simple_east_M7343ncDisplay.SetScaleArray = ['POINTS', 'z']
tanioka_simple_east_M7343ncDisplay.ScaleTransferFunction = 'PiecewiseFunction'
tanioka_simple_east_M7343ncDisplay.OpacityArray = ['POINTS', 'z']
tanioka_simple_east_M7343ncDisplay.OpacityTransferFunction = 'PiecewiseFunction'
tanioka_simple_east_M7343ncDisplay.DataAxesGrid = 'GridAxesRepresentation'
tanioka_simple_east_M7343ncDisplay.PolarAxes = 'PolarAxesRepresentation'
tanioka_simple_east_M7343ncDisplay.ScalarOpacityUnitDistance = 2050.9726164022213
tanioka_simple_east_M7343ncDisplay.ScalarOpacityFunction = zPWF
tanioka_simple_east_M7343ncDisplay.OpacityArrayName = ['POINTS', 'z']
tanioka_simple_east_M7343ncDisplay.SliceFunction = 'Plane'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tanioka_simple_east_M7343ncDisplay.ScaleTransferFunction.Points = [-0.8402088284492493, 0.0, 0.5, 0.0, 0.9494520425796509, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tanioka_simple_east_M7343ncDisplay.OpacityTransferFunction.Points = [-0.8402088284492493, 0.0, 0.5, 0.0, 0.9494520425796509, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
tanioka_simple_east_M7343ncDisplay.SliceFunction.Origin = [627500.0, 7339000.0, 0.0]

# show color bar/color legend
tanioka_simple_east_M7343ncDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(hFFtestsurfacexdmf)

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=hFFtestsurfacexdmf)
cellDatatoPointData1.CellDataArraytoprocess = ['U', 'V', 'W', 'partition', 'u', 'v', 'w']

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'
cellDatatoPointData1Display.ColorArrayName = ['POINTS', 'partition']
cellDatatoPointData1Display.LookupTable = partitionLUT
cellDatatoPointData1Display.SelectTCoordArray = 'None'
cellDatatoPointData1Display.SelectNormalArray = 'None'
cellDatatoPointData1Display.SelectTangentArray = 'None'
cellDatatoPointData1Display.OSPRayScaleArray = 'partition'
cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.SelectOrientationVectors = 'None'
cellDatatoPointData1Display.ScaleFactor = 20500.0
cellDatatoPointData1Display.SelectScaleArray = 'partition'
cellDatatoPointData1Display.GlyphType = 'Arrow'
cellDatatoPointData1Display.GlyphTableIndexArray = 'partition'
cellDatatoPointData1Display.GaussianRadius = 1025.0
cellDatatoPointData1Display.SetScaleArray = ['POINTS', 'partition']
cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.OpacityArray = ['POINTS', 'partition']
cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
cellDatatoPointData1Display.ScalarOpacityFunction = partitionPWF
cellDatatoPointData1Display.ScalarOpacityUnitDistance = 3530.290776620291
cellDatatoPointData1Display.OpacityArrayName = ['POINTS', 'partition']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
cellDatatoPointData1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
cellDatatoPointData1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(hFFtestsurfacexdmf, renderView1)

# show color bar/color legend
cellDatatoPointData1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=cellDatatoPointData1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.CoordinateResults = 1
calculator1.Function = 'coordsX*iHat + coordsY*jHat'

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'partition']
calculator1Display.LookupTable = partitionLUT
calculator1Display.SelectTCoordArray = 'None'
calculator1Display.SelectNormalArray = 'None'
calculator1Display.SelectTangentArray = 'None'
calculator1Display.OSPRayScaleArray = 'partition'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 20500.0
calculator1Display.SelectScaleArray = 'partition'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'partition'
calculator1Display.GaussianRadius = 1025.0
calculator1Display.SetScaleArray = ['POINTS', 'partition']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'partition']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.ScalarOpacityFunction = partitionPWF
calculator1Display.ScalarOpacityUnitDistance = 3530.1894499531313
calculator1Display.OpacityArrayName = ['POINTS', 'partition']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(cellDatatoPointData1, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Resample With Dataset'
resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=calculator1,
    DestinationMesh=tanioka_simple_east_M7343nc)
resampleWithDataset1.CellLocator = 'Static Cell Locator'

# show data in view
resampleWithDataset1Display = Show(resampleWithDataset1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
resampleWithDataset1Display.Representation = 'Slice'
resampleWithDataset1Display.ColorArrayName = ['POINTS', 'partition']
resampleWithDataset1Display.LookupTable = partitionLUT
resampleWithDataset1Display.SelectTCoordArray = 'None'
resampleWithDataset1Display.SelectNormalArray = 'None'
resampleWithDataset1Display.SelectTangentArray = 'None'
resampleWithDataset1Display.OSPRayScaleArray = 'partition'
resampleWithDataset1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleWithDataset1Display.SelectOrientationVectors = 'None'
resampleWithDataset1Display.ScaleFactor = 20500.0
resampleWithDataset1Display.SelectScaleArray = 'partition'
resampleWithDataset1Display.GlyphType = 'Arrow'
resampleWithDataset1Display.GlyphTableIndexArray = 'partition'
resampleWithDataset1Display.GaussianRadius = 1025.0
resampleWithDataset1Display.SetScaleArray = ['POINTS', 'partition']
resampleWithDataset1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.OpacityArray = ['POINTS', 'partition']
resampleWithDataset1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleWithDataset1Display.PolarAxes = 'PolarAxesRepresentation'
resampleWithDataset1Display.ScalarOpacityUnitDistance = 2050.9726164022213
resampleWithDataset1Display.ScalarOpacityFunction = partitionPWF
resampleWithDataset1Display.OpacityArrayName = ['POINTS', 'partition']
resampleWithDataset1Display.IsosurfaceValues = [6.5]
resampleWithDataset1Display.SliceFunction = 'Plane'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleWithDataset1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleWithDataset1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 13.0, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
resampleWithDataset1Display.SliceFunction.Origin = [627500.0, 7339000.0, 0.0]

# hide data in view
Hide(tanioka_simple_east_M7343nc, renderView1)

# hide data in view
Hide(calculator1, renderView1)

# show color bar/color legend
resampleWithDataset1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(tanioka_simple_east_M7343nc)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=tanioka_simple_east_M7343ncDisplay.SliceFunction)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=tanioka_simple_east_M7343ncDisplay)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=tanioka_simple_east_M7343ncDisplay.SliceFunction)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=tanioka_simple_east_M7343ncDisplay)

# set active source
SetActiveSource(resampleWithDataset1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=resampleWithDataset1Display.SliceFunction)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=resampleWithDataset1Display)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=resampleWithDataset1Display.SliceFunction)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=resampleWithDataset1Display)

# create a new 'Append Attributes'
appendAttributes1 = AppendAttributes(registrationName='AppendAttributes1', Input=[tanioka_simple_east_M7343nc, resampleWithDataset1])

# show data in view
appendAttributes1Display = Show(appendAttributes1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
appendAttributes1Display.Representation = 'Slice'
appendAttributes1Display.ColorArrayName = ['POINTS', 'z']
appendAttributes1Display.LookupTable = zLUT
appendAttributes1Display.SelectTCoordArray = 'None'
appendAttributes1Display.SelectNormalArray = 'None'
appendAttributes1Display.SelectTangentArray = 'None'
appendAttributes1Display.OSPRayScaleArray = 'U'
appendAttributes1Display.OSPRayScaleFunction = 'PiecewiseFunction'
appendAttributes1Display.SelectOrientationVectors = 'None'
appendAttributes1Display.ScaleFactor = 20500.0
appendAttributes1Display.SelectScaleArray = 'None'
appendAttributes1Display.GlyphType = 'Arrow'
appendAttributes1Display.GlyphTableIndexArray = 'None'
appendAttributes1Display.GaussianRadius = 1025.0
appendAttributes1Display.SetScaleArray = ['POINTS', 'U']
appendAttributes1Display.ScaleTransferFunction = 'PiecewiseFunction'
appendAttributes1Display.OpacityArray = ['POINTS', 'U']
appendAttributes1Display.OpacityTransferFunction = 'PiecewiseFunction'
appendAttributes1Display.DataAxesGrid = 'GridAxesRepresentation'
appendAttributes1Display.PolarAxes = 'PolarAxesRepresentation'
appendAttributes1Display.ScalarOpacityUnitDistance = 2050.9726164022213
appendAttributes1Display.ScalarOpacityFunction = zPWF
appendAttributes1Display.OpacityArrayName = ['POINTS', 'U']
appendAttributes1Display.SliceFunction = 'Plane'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
appendAttributes1Display.ScaleTransferFunction.Points = [-3.3200690514369273, 0.0, 0.5, 0.0, 2.813021656038931, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
appendAttributes1Display.OpacityTransferFunction.Points = [-3.3200690514369273, 0.0, 0.5, 0.0, 2.813021656038931, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
appendAttributes1Display.SliceFunction.Origin = [627500.0, 7339000.0, 0.0]

# hide data in view
Hide(resampleWithDataset1, renderView1)

# hide data in view
Hide(tanioka_simple_east_M7343nc, renderView1)

# show color bar/color legend
appendAttributes1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=appendAttributes1)
calculator2.Function = ''

# Properties modified on calculator2
calculator2.ResultArrayName = '$u_B$'
calculator2.Function = 'z-W'

# show data in view
calculator2Display = Show(calculator2, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'u_B'
u_BLUT = GetColorTransferFunction('u_B')

# get opacity transfer function/opacity map for 'u_B'
u_BPWF = GetOpacityTransferFunction('u_B')

# trace defaults for the display properties.
calculator2Display.Representation = 'Slice'
calculator2Display.ColorArrayName = ['POINTS', '$u_B$']
calculator2Display.LookupTable = u_BLUT
calculator2Display.SelectTCoordArray = 'None'
calculator2Display.SelectNormalArray = 'None'
calculator2Display.SelectTangentArray = 'None'
calculator2Display.OSPRayScaleArray = '$u_B$'
calculator2Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator2Display.SelectOrientationVectors = 'None'
calculator2Display.ScaleFactor = 20500.0
calculator2Display.SelectScaleArray = '$u_B$'
calculator2Display.GlyphType = 'Arrow'
calculator2Display.GlyphTableIndexArray = '$u_B$'
calculator2Display.GaussianRadius = 1025.0
calculator2Display.SetScaleArray = ['POINTS', '$u_B$']
calculator2Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator2Display.OpacityArray = ['POINTS', '$u_B$']
calculator2Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator2Display.DataAxesGrid = 'GridAxesRepresentation'
calculator2Display.PolarAxes = 'PolarAxesRepresentation'
calculator2Display.ScalarOpacityUnitDistance = 2050.9726164022213
calculator2Display.ScalarOpacityFunction = u_BPWF
calculator2Display.OpacityArrayName = ['POINTS', '$u_B$']
calculator2Display.IsosurfaceValues = [0.018896838735725507]
calculator2Display.SliceFunction = 'Plane'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator2Display.ScaleTransferFunction.Points = [-0.5217573833129869, 0.0, 0.5, 0.0, 0.5595510607844381, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator2Display.OpacityTransferFunction.Points = [-0.5217573833129869, 0.0, 0.5, 0.0, 0.5595510607844381, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
calculator2Display.SliceFunction.Origin = [627500.0, 7339000.0, 0.0]

# hide data in view
Hide(appendAttributes1, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1372, 763)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [627500.0, 7339000.0, 204177.80239484785]
renderView1.CameraFocalPoint = [627500.0, 7339000.0, 0.0]
renderView1.CameraViewAngle = 21.36304062909568
renderView1.CameraParallelScale = 137066.5896562689

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
