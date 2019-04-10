#!/usr/bin/env python2
#
# batch processing
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# January 2019

import MultiScaleData as MSD
import sys, yaml, os, shutil
from pytools import save_bin
from glob import glob
from itertools import product as comb, izip
import numpy as np
from chainmap import ChainMap
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser(
    description="Script producing the directory structure for ILLUMINA."
)
parser.add_argument("path", help="Path to the input parameters file [input_params.in].")
parser.add_argument("batch_name", nargs="?",
    help="Name for the produced batch file. "\
    "Will use the one defined in the parameter file if ommited.")

p = parser.parse_args()

def input_line(val,comment,n_space=30):
    value_str = ' '.join(str(v) for v in val)
    comment_str = ' ; '.join(comment)
    return "%-*s ! %s" % (n_space, value_str, comment_str )

def MSDOpen(filename,cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds

with open(os.path.join(p.path,"inputs_params.in")) as f:
    params = yaml.load(f,Loader=yaml.BaseLoader)

if p.batch_name is not None:
    params['batch_file_name'] = p.batch_name

for fname in glob(os.path.join(p.path,"%s*" % params['batch_file_name'])):
    os.remove(fname)

exp_name = params['exp_name']

# Add wavelength and multiscale
ds = MSDOpen("stable_lights.hdf5")

params['wavelength'] = np.loadtxt("wav.lst").tolist()
params['layer'] = range(len(ds))
params['observer_coordinates'] = zip(*ds.get_obs_pos())

for pname in ['layer','observer_coordinates']:
    if len(params[pname]) == 1:
        params[pname] = params[pname][0]

with open("zon.lst") as f:
    zones = f.read().split()

# Clear and create execution folder
dir_name = "exec" + os.sep
shutil.rmtree(dir_name,True)
os.makedirs(dir_name)

count = 0
multival = filter( lambda k: isinstance(params[k],list),params )
multival = sorted( multival, key=len, reverse=True ) # Semi-arbitrary sort
param_space = [ params[k] for k in multival ]
print "Number of executions:", np.prod(map(len,param_space))
for param_vals in comb(*param_space):
    local_params = OrderedDict(izip(multival,param_vals))
    P = ChainMap(local_params,params)
    if "azimuth_angle" in multival \
        and P['elevation_angle'] == 90 \
        and params['azimuth_angle'].index(P['azimuth_angle']) != 0:
        continue

    coords = P['observer_coordinates']
    if 'observer_coordinates' in multival:
        P['observer_coordinates'] = "%g_%g" % coords

    fold_name = dir_name + \
        os.sep.join(k+"_%s" % v for k,v in local_params.iteritems()) + os.sep

    os.makedirs(fold_name)
    #print fold_name

    wavelength = "%03d" % P["wavelength"]

    # Linking files
    mie_file = "%s_RH%02d_0.%s0um.mie.out" % (
        params['aerosol_profile'],
        params['relative_humidity'],
        wavelength )
    os.symlink(
        os.path.abspath(mie_file),
        fold_name+"aerosol.mie.out" )

    for zon in zones:
        os.symlink(
            os.path.abspath("fctem_wl_%s_zon_%s.dat" % (wavelength,zon)),
            fold_name+exp_name+"_fctem_"+zon+".dat" )

    ppath = os.environ['PATH'].split(os.pathsep)
    illumpath = filter(lambda s: "illumina" in s and "bin" in s, ppath)[0]
    os.symlink(
        os.path.abspath(illumpath+"/illumina"),
        fold_name+"illumina" )

    # Copying layer data
    layer = P["layer"]

    ds = MSDOpen("srtm.hdf5").extract_observer(coords)
    save_bin(fold_name+exp_name+"_topogra.bin", ds[layer])

    for name in ["obstd","obsth","obstf","altlp"]:
        ds = MSDOpen( "%s_%s.hdf5" % \
            ( exp_name, name ) ).extract_observer(coords)
        save_bin(fold_name+"%s_%s.bin" % \
            ( exp_name, name ), ds[layer] )

    for zon in zones:
        ds = MSDOpen( "%s_%s_lumlp_%s.hdf5" % \
            ( exp_name, wavelength, zon ) ).extract_observer(coords)
        ds.set_buffer(0)
        ds.set_overlap(0)
        save_bin(fold_name+"%s_lumlp_%s.bin" % \
            ( exp_name, zon ), ds[layer] )

    ds = MSDOpen( "modis_%s.hdf5" % \
        wavelength ).extract_observer(coords)
    save_bin(fold_name+"%s_reflect.bin" % \
        exp_name, ds[layer] )

    # Create illumina.in
    input_data = (
        (('', "Input file for ILLUMINA"),),
        ((exp_name, "Root file name"),),
        ((ds.pixel_size(layer), "Cell size along X [m]"),
         (ds.pixel_size(layer), "Cell size along Y [m]")),
        (("aerosol.mie.out", "Aerosol optical cross section file"),),
        (('', ''),),
        ((P['scattering_radius']*1000, "Double scattering radius [m]" ),
         (P['scattering_skip'], "Scattering step")),
        (('', ''),),
        ((wavelength, "Wavelength [nm]"),),
        ((P['air_pressure'], "Ground level pressure [kPa]"),),
        ((P['aerosol_optical_depth'], "500nm aerosol optical depth"),
         (P['angstrom_coefficient'], "Angstrom exponent")),
        ((len(zones), "Number of source types"),),
        ((P['stop_limit'], "Contribution threshold"),),
        (('', ''),),
        ((ds[layer].shape[1]/2 + 1, "Observer X position"),
         (ds[layer].shape[0]/2 + 1, "Observer Y position"),
         (P['observer_elevation'], "Observer elevation above ground [m]"),
         (1, "Beginning cell along line of sight")),
        (('', ''),),
        ((P['elevation_angle'], "Elevation viewing angle"),
         (P['azimuth_angle'], "Azimutal viewing angle")),
        (('', ''),),
        ((1., "Slit width [m]"),
         (1., "Pixel size [m]"),
         (1., "Focal length [m]"),
         (1.1283791671, "Apperture diameter [m]")),
        (('', "to get W/str/m**2, SW*PS*pi*AD**2/4/FL**2 = 1"),),
        (('', "1. 1. 1. 1.1283791671 is doing exactly that"),),
        ((P['cloud_model'], "Cloud model: "
            "0=clear, "
            "1=Thin Cirrus/Cirrostratus, "
            "2=Thick Cirrus/Cirrostratus, "
            "3=Altostratus/Altocumulus, "
            "4=Cumulus/Cumulonimbus, "
            "5=Stratocumulus"),),
        ((ds.pixel_size(layer)/2, "Minimal distance to nearest light source [m]"),)
    )

    with open(fold_name+"illumina.in",'w') as f:
        lines = ( input_line(*izip(*line_data)) for line_data in input_data )
        f.write( '\n'.join(lines) )

    # Write execute script
    with open(fold_name+"execute",'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --job-name=Illumina\n")
        f.write("#SBATCH --time=%d:00:00\n" % \
            params["estimated_computing_time"])
        f.write("#SBATCH --mem-per-cpu=1920\n")
        f.write("cd %s\n" % os.path.abspath(fold_name))
        f.write("umask 0011\n")
        f.write("./illumina\n")
    os.chmod(fold_name+"execute",0o777)

    # Append execution to batch list
    with open(
        sys.argv[1] + \
            '/' + \
            params['batch_file_name'] + \
            "_%d" % ((count/300)+1) ,
        'a' ) as f:
        f.write("cd %s\n" % os.path.abspath(fold_name))
        f.write("sbatch ./execute\n")
        f.write("sleep 0.05\n")

    count += 1

print "Final count:", count
print "Done."
