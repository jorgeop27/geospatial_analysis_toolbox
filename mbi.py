#!/usr/bin/env python

import gdal
import ia870
import mahotas as mh
import pymorph as pm
import numpy as np
from skimage import morphology,exposure,img_as_float,filters

def save_geotiff(outimage_file, image_array, base_geotiff, datatype = gdal.GDT_UInt16):
    rows = base_geotiff.RasterYSize
    cols = base_geotiff.RasterXSize
    geo_transform = base_geotiff.GetGeoTransform()
    geo_projection = base_geotiff.GetProjection()
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outimage_file, cols, rows, 1, datatype)
    outRaster.SetGeoTransform(geo_transform)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(image_array)
    outRaster.SetProjection(geo_projection)
    outband.FlushCache()
    outRaster = None
    return True

def opening_by_reconstruction(img, se):
    diffImg = True
    orig_img = img.copy()
    last_rec = img.copy()
    while diffImg:
        er_img = morphology.erosion(img, se)
        img = morphology.reconstruction(er_img, orig_img)
        if np.array_equal(last_rec, img):
            diffImg = False
        else:
            last_rec = img.copy()
    return last_rec

def proc_mbi(imgarray):
    # Normalize image:
    img = img_as_float(imgarray,force_copy=True)
    # Image equalization (Contrast stretching):
    p2,p98 = np.percentile(img, (2,98))
    img = exposure.rescale_intensity(img, in_range=(p2, p98), out_range=(0, 1))
    # Gamma correction:
    #img = exposure.adjust_gamma(img, 0.5)
    # Or Sigmoid correction:
    img = exposure.adjust_sigmoid(img)
    
    print "Init Morph Proc..."
    sizes = range(2,40,5)
    angles = [0,18,36,54,72,90,108,126,144,162]
    szimg = img.shape
    all_thr = np.zeros((len(sizes),szimg[0], szimg[1])).astype('float64')
    all_dmp = np.zeros((len(sizes) - 1,szimg[0], szimg[1])).astype('float64')
    
    idx = 0
    for sz in sizes:
        print sz
        builds_by_size = np.zeros(szimg).astype('float64')
        for ang in angles:
            print ang
            stel = ia870.iaseline(sz, ang)
            oprec = opening_by_reconstruction(img, stel)
            thr = np.absolute(img-oprec)
            builds_by_size += thr
        all_thr[idx,:,:] = (builds_by_size / len(angles))
        if idx>0:
            all_dmp[idx-1,:,:] = all_thr[idx,:,:] - all_thr[idx-1,:,:]
        idx += 1
    mbi = np.mean(all_dmp, axis=0)
    return mbi

def main(filename,outfile):
    # Open geotiff file:
    geoimg = gdal.Open(filename, gdal.GA_ReadOnly)
    image = geoimg.ReadAsArray().astype(np.uint16)
    rows = geoimg.RasterYSize
    cols = geoimg.RasterXSize
    if geoimg.RasterCount>1:
        img = image.max(axis=0)
    else:
        img = image.copy();
    mbi = proc_mbi(img)
    save_geotiff(outfile, mbi, geoimg, gdal.GDT_Float64)
    # otsu = filters.threshold_otsu(mbi,1024)
    # mbi_bw = (mbi > otsu)
    # save_geotiff('/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbiBW_otsu.tif', mbi_bw, geoimg, gdal.GDT_Byte)
    # yen = filters.threshold_yen(mbi,1024)
    # mbi_bw = (mbi > yen)
    # save_geotiff('/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbiBW_yen.tif', mbi_bw, geoimg, gdal.GDT_Byte)
    # isodata = filters.threshold_isodata(mbi,1024)
    # mbi_bw = (mbi > isodata)
    # save_geotiff('/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbiBW_isodata.tif', mbi_bw, geoimg, gdal.GDT_Byte)
    li = filters.threshold_li(mbi)
    mbi_bw = (mbi > li)
    save_geotiff('/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbiBW_li.tif', mbi_bw, geoimg, gdal.GDT_Byte)
    # mbi_bw = (mbi > 0)
    # save_geotiff('/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbiBW_minVal.tif', mbi_bw, geoimg, gdal.GDT_Byte)
    return True

if __name__ == '__main__':
    ## Main Function:
    # image = 'Z:/Documents/LSE/Sample1m/Sample1m.tif'
    # shapefile = 'Z:/Documents/LSE/Sample1m/Sample1m_Roads.shp'
    # dmpfile = 'Z:/Documents/LSE/Sample1m/Sample1m_dmp_output/Sample1m_mbiBW.tif'
    # cluster_file = 'Z:/Documents/LSE/Sample1m/Sample1m_Roads_raster_cluster.tif'
    image = '/Users/jorge/Documents/LSE/Sample1m/Sample1m.tif'
    output = '/Users/jorge/Documents/LSE/Sample1m/Sample1m_mbi.tif'
    # shapefile = '/Users/jorge/Documents/LSE/Sample_1m/Sample_Pan1m_Roads.shp'
    main(image, output)
