try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

combxyz_bin = PLOT3DReader( FileName='/data/turbine_Stg/zDIR.P3D.rel.6201-11001/s35_noinj.r2b1.p3d.g6201' )
combxyz_bin.Functions = [110]
combxyz_bin.QFileName = ['/data/turbine_Stg/zDIR.P3D.rel.6201-11001/s35_noinj.r2b1.p3d.q6201']
combxyz_bin.FunctionFileName = ''

Render()
