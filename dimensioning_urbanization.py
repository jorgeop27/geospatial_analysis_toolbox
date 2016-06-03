import arcpy
from arcpy import env
import numpy as np
import math

imgDir = "Z:\Documents\LSE\\2013\\"
env.workspace = imgDir
env.overwriteOutput = True

def main(file, radio_length=5):
    # Temporary features
    binary_shapefile = "temp_binary.shp"
    distances_dbf = "temp_distances.dbf"
    # Outputs
    centroids = "nodes2.shp"
    edges = "edges2.shp"
    
    # Convert binary raster to polygon and calculate geometry attributes
    arcpy.RasterToPolygon_conversion(file, binary_shapefile)
    arcpy.AddGeometryAttributes_management(binary_shapefile, "AREA;PERIMETER_LENGTH;CENTROID", "METERS", "SQUARE_METERS")

    # Calculate tree nodes (polygon centroids)
    arcpy.FeatureToPoint_management(binary_shapefile, centroids, "CENTROID")
    # Calculate distances with a radial length of 'radio_length':
    arcpy.PointDistance_analysis(centroids, centroids, distances_dbf, str(radio_length))
    
    # Create a new shapefile for the edges of the graph
    arcpy.CreateFeatureclass_management(env.workspace, edges, "POLYLINE", spatial_reference=centroids)
    
    # Add fields to calculate to nodes and edges shapefiles
    arcpy.AddField_management(centroids, "SHAPE_IND", "FLOAT")
    arcpy.AddField_management(centroids, "DEGREE", "SHORT")
    arcpy.AddField_management(centroids, "NEARNESS", "FLOAT")
    arcpy.AddField_management(edges, "EDGE_LEN", "FLOAT")
    arcpy.AddField_management(edges, "SIGNIFICAN", "FLOAT")
    
    # Calculate shape index
    arcpy.CalculateField_management(centroids, "SHAPE_IND", "!PERIMETER! / (4*math.sqrt(!POLY_AREA!))", "PYTHON_9.3")

    # Create cursor to add new polylines (Edges of the graph) and it's attributes (Length and Local significance)
    cursor = arcpy.da.InsertCursor(edges, ["SHAPE@","EDGE_LEN","SIGNIFICAN"])

    # Loop through the polygons in centroids shapefile and update each centroid with its calculated properties
    points = arcpy.Array()
    # rows = arcpy.da.UpdateCursor(centroids,["ORIG_FID","CENTROID_X","CENTROID_Y","POLY_AREA","DEGREE","NEARNESS"])
    rows = arcpy.da.SearchCursor(distances_dbf, ["INPUT_FID","NEAR_FID","DISTANCE"])

    dists = {}
    for rw in rows:
        # Find the centroid point, the original point of the edge
        objId = rw[0]
        ngId = rw[1]
        length = rw[2]
        where = "ORIG_FID = %d" % objId
        node_cursor = arcpy.da.UpdateCursor(centroids,["CENTROID_X","CENTROID_Y","POLY_AREA","DEGREE"],where)
        node = node_cursor.next()
        org_point = arcpy.Point()
        org_point.X = node[0]
        org_point.Y = node[1]
        points.add(org_point)
        # Find the area of the centroid's polygon
        org_area = node[2]
        node[3] = (node[3] + 1)
        if objId not in dists.keys():
            dists[objId] = [length]
        else:
            dists[objId].append(length)
        node_cursor.updateRow(node)
        del node_cursor
        # Find the neighbour nodes
        where = "ORIG_FID = %d" % ngId
        ngNode_cursor = arcpy.da.SearchCursor(centroids,["CENTROID_X","CENTROID_Y","POLY_AREA"],where)
        ngNode = ngNode_cursor.next()
        # Create the destiny point of the edge, add it to the points array and create the polyline
        ng_point = arcpy.Point()
        ng_point.X = ngNode[0]
        ng_point.Y = ngNode[1]
        points.add(ng_point)
        line = arcpy.Polyline(points)
        # Find the area of the neighbour
        ng_area = ngNode[2]
        del ngNode_cursor
        # Calculate local significance
        significance = (org_area*ng_area)/math.pow(length,2)
        # Insert new polyline (the edge) and calculated attributes
        cursor.insertRow([line,length,significance])
        # Clearn the points array for the next loop
        points.removeAll()
    del cursor, rows
    for objId in dists.keys():
        where = "ORIG_FID = %d" % objId
        node_cursor = arcpy.da.UpdateCursor(centroids,["NEARNESS"],where)
        node = node_cursor.next()
        near_distances = np.array(dists[objId])
        # Calculate local nearness
        near_distances.sort()
        nearness = np.sum(near_distances[0:3])
        # Update the node with the calculated attributes
        node[0] = nearness
        node_cursor.updateRow(node)
        del node_cursor
    where = "DEGREE = 0"
    cursor = arcpy.da.UpdateCursor(centroids,["DEGREE"],where)
    for row in cursor:
        cursor.deleteRow()
    del cursor
    

### MAIN ###
main('PAN_mbiBW.tif', radio_length=10)