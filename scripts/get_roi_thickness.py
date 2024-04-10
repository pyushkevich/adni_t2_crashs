#!/usr/bin/env python3
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from statsmodels.stats.weightstats import DescrStatsW
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# ROI names
roi_names = {
    0: 'White_Matter',
    1: 'ERC',
    2: 'BA35',
    3: 'BA36',
    4: 'PHC'
}

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('subject', help='Subject ID')
parser.add_argument('session', help='Session ID')
parser.add_argument('scan', help='Scan ID')
parser.add_argument('side', help='Side (left/right)')
parser.add_argument('mesh', help='Mesh from which to read thickness')
parser.add_argument('output', help='Output file with the statistics')
args = parser.parse_args()

# Read the mesh
r = vtk.vtkPolyDataReader()
r.SetFileName(args.mesh)
r.Update()
mesh = r.GetOutput()
tri=vtk_to_numpy(mesh.GetPolys().GetData()).reshape(-1,4)[:,1:]
x=vtk_to_numpy(mesh.GetPoints().GetData())

# Compute triangle areas
u = x[tri[:,1],:]-x[tri[:,0],:]
v = x[tri[:,2],:]-x[tri[:,0],:]
area = np.sqrt(np.sum(np.cross(u,v)**2,1)) / 2

# Read thickness and label arrays
label = vtk_to_numpy(mesh.GetCellData().GetArray('label'))
label_ids = np.unique(label)

# Get the triangle area for each label
tri_area_by_label = np.zeros((label.shape[0], len(label_ids)))
for i,k in enumerate(label_ids):
    tri_area_by_label[label==k, i] = area[label==k]

# Get the triangle area for each vertex
vtx_area_by_label = np.zeros((x.shape[0], len(label_ids)))
for j in range(3):
    vtx_area_by_label[tri[:,j]] += tri_area_by_label / 3

# Compute the weighted statistics for the cm-rep based thickness computation
thick_voronoi = vtk_to_numpy(mesh.GetPointData().GetArray('VoronoiRadius'))
m_voronoi = ~np.isnan(thick_voronoi)
stat_voronoi = { k: DescrStatsW(thick_voronoi[m_voronoi], vtx_area_by_label[m_voronoi,i]) for (i,k) in enumerate(label_ids) }

# Compute the weighted statistics for the Cruise-based thickness computation
thick_cruise_times_mask = vtk_to_numpy(mesh.GetPointData().GetArray('CruiseThicknessTimesMask'))
thick_cruise_mask = vtk_to_numpy(mesh.GetPointData().GetArray('CruiseThicknessMask'))
m_cruise = thick_cruise_mask >= 0.5
thick_cruise = thick_cruise_times_mask[m_cruise] / thick_cruise_mask[m_cruise]
stat_cruise = { k: DescrStatsW(thick_cruise, vtx_area_by_label[m_cruise,i]) for (i,k) in enumerate(label_ids) }

# Compute the dataframe with the statistics
df = pd.DataFrame({
    'subject': args.subject,
    'session': args.session,
    'scan': args.scan,
    'side': args.side,
    'label': [ roi_names[k] for k in label_ids ],
    'VR_mean': [ stat_voronoi[k].mean for k in label_ids ],
    'VR_median': [ stat_voronoi[k].quantile(0.5, False)[0] for k in label_ids ],
    'VR_q95': [ stat_voronoi[k].quantile(0.95, False)[0] for k in label_ids ],
    'CT_mean': [ stat_cruise[k].mean for k in label_ids ],
    'CT_median': [ stat_cruise[k].quantile(0.5, False)[0] for k in label_ids ],
    'CT_q95': [ stat_cruise[k].quantile(0.95, False)[0] for k in label_ids ]
})

# Save the statistics to the output file
df.to_csv(args.output)

