#!MC 1410
$!SetContourVar 
  Var = 14
  ContourGroup = 1
  LevelInitMode = ResetToNice
$!GlobalContour 1  ColorMapName = 'Small Rainbow'
$!ContourLevels New
  ContourGroup = 1
  RawData
17
0.01
0.028125
0.04625
0.064375
0.0825
0.100625
0.11875
0.136875
0.155
0.173125
0.19125
0.209375
0.2275
0.245625
0.26375
0.281875
0.3
$!SetContourVar 
  Var = 4
  ContourGroup = 2
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 5
  ContourGroup = 3
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 6
  ContourGroup = 4
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 7
  ContourGroup = 5
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 8
  ContourGroup = 6
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 9
  ContourGroup = 7
  LevelInitMode = ResetToNice
$!SetContourVar 
  Var = 10
  ContourGroup = 8
  LevelInitMode = ResetToNice
$!IsoSurfaceLayers Show = Yes
$!IsoSurfaceAttributes 1  Isovalue1 = 0.0700000000000000067
$!RedrawAll 
$!Pick DeselectAll
$!Pick AddAllInRect
  SelectText = Yes
  SelectGeoms = Yes
  SelectZones = Yes
  ConsiderStyle = Yes
  X1 = 2.68172983479
  X2 = 7.48639455782
  Y1 = 3.19655004859
  Y2 = 7.28595724004
