import numpy as np
import arcpy
import os
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt
# Imports for Spider chart:
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection


imgDir = "Z:\Documents\LSE"
arcpy.env.workspace = imgDir
arcpy.env.overwriteOutput = True
if not arcpy.CheckOutExtension('Spatial'):
    print "This application needs the Spatial Analysis Toolbox"


# Functions for Spider (Radar) Chart:
def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.
    This function creates a RadarAxes projection and registers it.
    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.
    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
    # rotate theta such that the first axis is at the top
    theta += np.pi/2
    
    def draw_poly_patch(self):
        verts = unit_poly_verts(theta)
        return plt.Polygon(verts, closed=True, edgecolor='k')
    
    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)
    
    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)
    
    class RadarAxes(PolarAxes):
        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]
        
        def fill(self, *args, **kwargs):
            """Override fill so that line is closed by default"""
            closed = kwargs.pop('closed', True)
            return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)
        
        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super(RadarAxes, self).plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)
        
        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)
        
        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)
        
        def _gen_axes_patch(self):
            return self.draw_patch()
        
        def _gen_axes_spines(self):
            if frame == 'circle':
                return PolarAxes._gen_axes_spines(self)
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.
            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)
            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}
    register_projection(RadarAxes)
    return theta

def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.
    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

# End of Spider Chart functions
##############################################################################

def plot_charts(files, output_files):
    cmap = plt.cm.get_cmap('Spectral')
    colors = [cmap(i) for i in np.linspace(0,1,(len(output_files)+1))]
    labels = [os.path.split(os.path.splitext(z)[0])[1] for z in files]
    full_data = []
    for of in output_files:
        data = np.array([x for x in arcpy.da.SearchCursor(of, ["BUILD_DENS","AVG_SIZE","HETEROGENE",
                                                            "AVG_DIST1","AVG_DIST2","AVG_DIST3","AVG_DIST4",
                                                            "AVG_DIST5","AVG_DIST6","AVG_DIST7","AVG_DIST8",
                                                            "AVG_DIST9","AVG_DIST10"])])
        full_data.append(data)
    # BoxPlot of Building density:
    f1 = plt.figure(1)
    ax1 = f1.add_subplot(111)
    density = []
    for d in full_data:
        density.append(d[:,0])
    bp1 = ax1.boxplot(density, patch_artist=True)
    for col_idx,box in enumerate(bp1['boxes']):
        box.set(facecolor = colors[col_idx+1], linewidth=2)
    for flier in bp1['fliers']:
        flier.set(marker='o', color='#ff298a', alpha=0.7)
    ax1.set_xticklabels(labels)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    plt.title('Box Plot of Building Density')
    plt.ylabel('Building Density')
    # BoxPlot of Average Buildings sizes:
    f2 = plt.figure(2)
    ax2 = f2.add_subplot(111)
    avg_size = []
    for d in full_data:
        avg_size.append(d[:,1])
    bp2 = ax2.boxplot(avg_size, patch_artist=True)
    for col_idx,box in enumerate(bp2['boxes']):
        box.set(facecolor = colors[col_idx+1], linewidth=2)
    for flier in bp2['fliers']:
        flier.set(marker='o', color='#ff298a', alpha=0.7)
    ax2.set_xticklabels(labels)
    ax2.get_xaxis().tick_bottom()
    ax2.get_yaxis().tick_left()
    plt.title('Box Plot of Avg. Building Size')
    plt.ylabel('Average Building Size')
    # # BoxPlot of Heterogeneity:
    f3 = plt.figure(3)
    ax3 = f3.add_subplot(111)
    heterogeneity = []
    for d in full_data:
        heterogeneity.append(d[:,2])
    bp3 = ax3.boxplot(heterogeneity, patch_artist=True)
    for col_idx,box in enumerate(bp3['boxes']):
        box.set(facecolor = colors[col_idx+1], linewidth=2)
    for flier in bp3['fliers']:
        flier.set(marker='o', color='#ff298a', alpha=0.7)
    ax3.set_xticklabels(labels)
    ax3.get_xaxis().tick_bottom()
    ax3.get_yaxis().tick_left()
    plt.title('Box Plot of Heterogeneity')
    plt.ylabel('Heterogeneity')
    # Plot of Average Building distances:
    f4 = plt.figure(4)
    # avg_distances = []
    for idx,d in enumerate(full_data):
        avg_distances = d[:,3:].mean(axis=0)
        plt.plot(range(1,11), avg_distances, color=colors[idx+1], linestyle='-', linewidth=2, marker='^', \
                 markersize=10, markerfacecolor=colors[idx+1], label=labels[idx])
    plt.xlim(0,11)
    plt.yscale('log')
    plt.grid(True)
    plt.ylabel('Average Building Distance (Log scale)')
    plt.xlabel('Number of Building Neighbours')
    plt.legend()
    # Spider Chart of Density, Size, Heterogeneity and Distance to 1 nearest neighbour
    N = 4
    theta = radar_factory(N)#, frame='polygon')
    spoke_labels = ['BD','ABS','HI','ABD']
    f5 = plt.figure()
    ax5 = f5.add_subplot(111, projection='radar')
    for idx,d in enumerate(full_data):
        data_radar = d[:,0:4]
        data_radar = (data_radar / data_radar.max(axis=0)).mean(axis=0)
        ax5.plot(theta, data_radar, color=colors[idx+1])
        ax5.fill(theta, data_radar, facecolor=colors[idx+1], alpha=0.25)
    ax5.set_varlabels(spoke_labels)
    plt.legend(labels, loc=(0.65, 0.94), labelspacing=0.1)
    plt.show()

def indices_calculation(buildings, blocks_file=False, block_size=300):
    # This function calculates the Building Density, Average Building Size,
    # Average distance to nearest neighbour (1 up to 10 neighbours) and
    # the Heterogeneity index from a buildings shape file
    #
    # - buildings is the file name of the buildings. It should be a shapefile of polygons
    # - blocks_file is an optional shapefile of district boundaries or grids. This file will
    # will be used for dividing the buildings file and compute the statistics for each block.
    # The results will also be put in this file's data table. If this parameter is not input,
    # the script creates automatically a grid file of block_size x block_size per grid
    # - block_size is the size of the grid to use in case no blocks_file is input.
    #
    
    temp_dir = "__temp_output"
    if os.path.exists(os.path.join(imgDir,temp_dir)):
        shutil.rmtree(os.path.join(imgDir,temp_dir))
    os.makedirs(os.path.join(imgDir,temp_dir))
    desc = arcpy.Describe(buildings)
    extents = desc.extent

    net_orig = "%f %f" % (extents.XMin, extents.YMin)
    net_oppos = "%f %f" % (extents.XMax, extents.YMax)
    net_orient = "%f %f" % (extents.XMin, extents.YMin+1)
    if blocks_file==False:
        basename,ext = os.path.splitext(buildings)
        blocks_file = '%s_output.shp' % basename
        arcpy.CreateFishnet_management(blocks_file, net_orig, net_orient, block_size, block_size, 0, 0, net_oppos, "NO_LABELS", "#", "POLYGON")
    
    # Add a new field to the blocks_file for use in the Split_analysis function:
    arcpy.AddField_management(blocks_file, "FID_TEXT", "TEXT")
    # Cast FID as Text and save it in the new field:
    arcpy.CalculateField_management(blocks_file, "FID_TEXT", "str(!FID!)", "PYTHON_9.3")
    # Find blocks neighbours:
    neighbours_file = os.path.join(temp_dir,"blocks_neighbours.dbf")
    arcpy.PolygonNeighbors_analysis(blocks_file, neighbours_file, ["FID"], True)
    # Split original buildings file using the blocks file:
    arcpy.Split_analysis(buildings, blocks_file, "FID_TEXT", temp_dir)
    
    # Add new fields for calculated indices:
    arcpy.AddField_management(blocks_file, "BUILD_DENS", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_SIZE", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST1", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST2", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST3", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST4", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST5", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST6", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST7", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST8", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST9", "FLOAT")
    arcpy.AddField_management(blocks_file, "AVG_DIST10", "FLOAT")
    arcpy.AddField_management(blocks_file, "HETEROGENE", "FLOAT")
    
    with arcpy.da.UpdateCursor(blocks_file,["FID_TEXT","SHAPE@AREA","BUILD_DENS","AVG_SIZE",
                                            "AVG_DIST1","AVG_DIST2","AVG_DIST3","AVG_DIST4",
                                            "AVG_DIST5","AVG_DIST6","AVG_DIST7","AVG_DIST8",
                                            "AVG_DIST9","AVG_DIST10"]) as rows:
        for rw in rows:
            mem_grid = os.path.join(temp_dir,"%s.shp" % rw[0])
            if not os.path.exists(os.path.join(imgDir,mem_grid)):
                continue
            stats_file = os.path.join(temp_dir,"%s_stats.dbf" % rw[0])
            arcpy.Statistics_analysis(mem_grid, stats_file, [["Shape_Area","SUM"],["Shape_Area","MEAN"]])
            # Fields: [FREQUENCY, SUM_Shape_, MEAN_Shape]
            with arcpy.da.SearchCursor(stats_file, ["SUM_Shape_","MEAN_Shape"]) as line:
                line_data = line.next()
            rw[2] = line_data[0]/rw[1]
            rw[3] = line_data[1]
            # Get the centroids of the polygons to calculate the distances to neighbours:
            centroids = [x[0] for x in arcpy.da.SearchCursor(mem_grid, "SHAPE@XY")]
            # Convert centroids to a numpy array of complex numbers and calculate array of distances:
            cents = np.array([[complex(c[0],c[1]) for c in centroids]])
            distance_array = abs(cents.T - cents)
            # Sort each row of the array (increasing)
            distance_array.sort(axis=1)
            # Calculate average distance for n = 1:10 neighbours:
            for n in xrange(10):
                avg_dist = distance_array[:,1:(2+n)].mean()
                if np.isnan(avg_dist):
                    avg_dist = 0
                rw[4+n] = avg_dist
            rows.updateRow(rw)
    
    # Calculate Heterogeneity:
    densities = {fid:density for (fid,density) in arcpy.da.SearchCursor(blocks_file, ["FID", "BUILD_DENS"])}
    with arcpy.da.UpdateCursor(blocks_file,["FID","HETEROGENE"]) as cursor:
        for rw in cursor:
            where = "src_FID=%d" % rw[0]
            block_density = densities[rw[0]]
            neighbours = [x[0] for x in arcpy.da.SearchCursor(neighbours_file,["nbr_FID"],where)]
            sum_density = 0.0
            for nbr in neighbours:
                sum_density += abs(block_density - densities[nbr])
            heterogeneity = sum_density/len(neighbours)
            rw[1] = heterogeneity
            cursor.updateRow(rw)
    return blocks_file

# files is an array of buildings files (shapefiles of polygons)
files = ['TrueBuildings_Sample.shp','TrueBuildings_Sample2.shp','TrueBuildings_Sample3.shp']
# block_files is an array of blocks (or district boundaries) to use for the calculations.
# Each file in the 'files' array must have it's corresponding block file in the 'block_files' array.
# If block_files is empty, the script creates a grid of 300x300 for each file.
block_files = []
# output_files will store the names of the files with the statistics generated
output_files = []
for idx, file in enumerate(files):
    if len(block_files)>0:
        output_file = indices_calculation(file, block_files[idx])
    else:
        output_file = indices_calculation(file, block_size=300)
    output_files.append(output_file)
# plot_charts generates the charts using the statistics calculated before
plot_charts(files,output_files)
