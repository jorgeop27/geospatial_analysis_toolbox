import numpy as np
import arcpy
import os
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

def plot_charts(file):
    data = np.array([x for x in arcpy.da.SearchCursor(file, ["BUILD_DENS","AVG_SIZE","HETEROGENE",
                                                        "AVG_DIST1","AVG_DIST2","AVG_DIST3","AVG_DIST4",
                                                        "AVG_DIST5","AVG_DIST6","AVG_DIST7","AVG_DIST8",
                                                        "AVG_DIST9","AVG_DIST10"])])
    # BoxPlot of Building density:
    density = data[:,0]
    plt.boxplot(density)
    plt.title('Box Plot of Building Density')
    plt.ylabel('Building Density')
    # BoxPlot of Average Buildings sizes:
    plt.figure()
    avg_size = data[:,1]
    plt.boxplot(avg_size)
    plt.title('Box Plot of Avg. Building Size')
    plt.ylabel('Average Building Size')
    # BoxPlot of Heterogeneity:
    plt.figure()
    heterogeneity = data[:,2]
    plt.boxplot(heterogeneity)
    plt.title('Box Plot of Heterogeneity')
    plt.ylabel('Heterogeneity')
    # Plot of Average Building distances:
    plt.figure()
    avg_distances = data[:,3:].mean(axis=0)
    plt.plot(range(1,11), avg_distances,'b^',range(1,11),avg_distances,'b')
    plt.xlim(0,11)
    plt.yscale('log')
    # plt.axis([0, 11, -10, 100])
    plt.grid(True)
    plt.ylabel('Average Building Distance (Log scale)')
    plt.xlabel('Number of Building Neighbours')
    # Spider Chart of Density, Size, Heterogeneity and Distance to 1 nearest neighbour
    data_radar = data[:,0:4]
    # Normalize data and calculate mean
    data_radar = (data_radar / data_radar.max(axis=0)).mean(axis=0)
    # print data_radar
    N = 4
    theta = radar_factory(N, frame='polygon')
    spoke_labels = ['BD','ABS','HI','ABD']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='radar')
    ax.plot(theta, data_radar, color='b')
    ax.fill(theta, data_radar, facecolor='b', alpha=0.25)
    ax.set_varlabels(spoke_labels)
    plt.show()

def indices_calculation(buildings, blocks_file=False, block_size=300):
    temp_dir = "temp_output"
    if not os.path.exists(os.path.join(imgDir,temp_dir)):
        os.makedirs(os.path.join(imgDir,temp_dir))
    desc = arcpy.Describe(buildings)
    extents = desc.extent

    net_orig = "%f %f" % (extents.XMin, extents.YMin)
    net_oppos = "%f %f" % (extents.XMax, extents.YMax)
    net_orient = "%f %f" % (extents.XMin, extents.YMin+1)
    if blocks_file==False:
        blocks_file = 'output_grid.shp'
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
                rw[4+n] = avg_dist
            rows.updateRow(rw)
    
    # Calculate Heterogeneity:
    densities = {fid:density for (fid,density) in arcpy.da.UpdateCursor(blocks_file, ["FID", "BUILD_DENS"])}
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


file = 'TrueBuildings_Sample.shp'
output_file = indices_calculation(file, block_size=200)
plot_charts(output_file)