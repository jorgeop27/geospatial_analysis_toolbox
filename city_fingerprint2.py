import arcpy
from arcpy import env
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage,dendrogram
import math
import os.path
import csv

shpDir = "Z:\\Documents\\LSE\\Metropolitan\\"
env.workspace = shpDir
env.overwriteOutput = True
if not arcpy.CheckOutExtension('Spatial'):
    print "This application needs the Spatial Analysis Toolbox"

def fingerprint(file, dphi, city_index=0, plot=False):
    [basename, ext] = os.path.splitext(file)
    outFile = "%s_blocks.shp" % basename
    boundFile = '%s_blocks_bounds.shp' % basename
    csvblocks = os.path.join(shpDir,'%s_blocks_metrics.csv' % basename)
    csvdistrib = os.path.join(shpDir,'%s_formfactor_distribution.csv' % basename)
    
    arcpy.FeatureToPolygon_management([file], outFile)
    arcpy.CalculateField_management(outFile, "Id", "!FID!", "PYTHON_9.3")
    arcpy.AddGeometryAttributes_management(outFile, ["AREA_GEODESIC","PERIMETER_LENGTH_GEODESIC"])
    arcpy.AddField_management(outFile, "SHAPE_IND", "FLOAT")
    arcpy.CalculateField_management(outFile, "SHAPE_IND", "!PERIM_GEO! / (4*math.sqrt(abs(!AREA_GEO!)))", "PYTHON_9.3")
    arcpy.MinimumBoundingGeometry_management(outFile, boundFile, 'CIRCLE', mbg_fields_option=True)
    arcpy.AddField_management(boundFile, "ORIG_AREA", "FLOAT")
    arcpy.CalculateField_management(boundFile, "ORIG_AREA", "!AREA_GEO!", "PYTHON_9.3")
    arcpy.AddGeometryAttributes_management(boundFile, "AREA_GEODESIC")
    arcpy.AddField_management(boundFile, "FORM_FACT", "FLOAT")
    arcpy.CalculateField_management(boundFile, "FORM_FACT", "!ORIG_AREA!/!AREA_GEO!", "PYTHON_9.3")
    
    with open(csvblocks, 'wb') as f1:
        wrt1 = csv.writer(f1)
        data = arcpy.da.TableToNumPyArray(boundFile, ["ORIG_FID","ORIG_AREA","PERIM_GEO","FORM_FACT","SHAPE_IND"])
        header1 = ['FID', 'AREA', 'PERIMETER', 'FORM_FACTOR', 'SHAPE_INDEX']
        wrt1.writerow(header1)
        wrt1.writerows(data)
    if plot:
        fig = plt.figure(city_index)
    area = data["ORIG_AREA"]
    area[area<0] = 0
    ffact = data['FORM_FACT']
    ffact[ffact<0] = 0
    ffact[ffact>=1] = 1
    total_cell_number = ffact.size
    with open(csvdistrib, 'wb') as f2:
        wrt2 = csv.writer(f2)
        header2 = ['AREA_TYPE']
        bins = np.linspace(0,1,dphi+1)
        for b in xrange(1,len(bins)):
            header2.append("%.2f<Phi<%.2f"%(bins[b-1],bins[b]))
        wrt2.writerow(header2)
        if plot:
            [n, bins, patches] = plt.hist(ffact, bins=bins, histtype='step', alpha=0.5, color='b', label='All Areas')
        else:
            [n, bins] = np.histogram(ffact, bins=bins)
        rw = ['ALL']
        rw.extend(list(n))
        wrt2.writerow(rw)
        nbins = 4
        org = [0, 1000, 10000, 100000]
        end = [1000, 10000, 100000, area.max()]
        colors = ['r','g','c','m','y','k']
        
        for idx in xrange(nbins):
            if idx==0:
                ff_indices = (area<end[idx])
                lbl = 'AREA<%d' % (end[idx],)
            elif idx==(nbins-1):
                lbl = 'AREA>%d' % (org[idx],)
                ff_indices = (area>=org[idx])
            else:
                ff_indices = (area>=org[idx])&(area<end[idx])
                lbl = '%d<AREA<%d' % (org[idx],end[idx])
            ffact_area = ffact[ff_indices]
            if ffact_area.size==0:
                continue
            alpha, bins = np.histogram(ffact_area, bins=bins)
            rw = [lbl]
            rw.extend(list(alpha))
            wrt2.writerow(rw)
            if idx==1:
                # Area between 10^3 and 10^4
                alpha1 = alpha.astype('double')/total_cell_number
                # print alpha1
            if idx==2:
                # Area between 10^4 and 10^5
                alpha2 = alpha.astype('double')/total_cell_number
                # print alpha2
            if plot:
                plt.hist(ffact_area, bins=bins, histtype='step', color=colors[idx], label=lbl)
    if plot:
        plt.legend(loc='upper right')
        plt.title('Distribution of Shape Factor\n%s'%file)
    return [alpha1, alpha2]

def city_comparison(file_array):
    dphi = 40
    step = 1.0/dphi
    cities_fingerprint = np.zeros([len(file_array),2,dphi],'double')
    indx = 0
    for file in file_array:
        [alph1,alph2] = fingerprint(file, dphi, indx, plot=False)
        # print alph1,alph2
        cities_fingerprint[indx,0,:] = alph1
        cities_fingerprint[indx,1,:] = alph2
        indx+=1
    distance_matrix = np.zeros([len(file_array),len(file_array)],'double')
    for i in xrange(len(file_array)):
        city_i = cities_fingerprint[i,:,:]
        d_cities = np.square(cities_fingerprint - city_i)*0.01
        d_cities = np.sum(d_cities,2)
        D_cities = np.sum(np.square(d_cities),1)
        distance_matrix[i,:] = D_cities
    distance_condensed = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(distance_condensed)
    g = plt.figure(indx)
    dendrogram(linkage_matrix, orientation='right')
    plt.title('Dendrogram\nHierarchical Clustering of %d Cities'%len(file_array))
    plt.show()
    return True

### MAIN ###
city_comparison(['SmallTest1.shp','SmallTest2.shp','SmallTest3.shp','SmallTest4.shp','SmallTest5.shp','SmallTest6.shp'])
# city_comparison(['johannesburg/johannesburg_south-africa_osm_roads.shp',''])