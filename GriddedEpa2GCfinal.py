#!/usr/bin/env python
# coding: utf-8

# # Preparation of GEOS-Chem Emissions from CMAQ
# 
#     author: Barron H. Henderson
#     contributors: Lyssa Freese. Adding names pending approval...
#     date: 2020-06-04

# # Overview
# 
# This notebook takes SMOKE outputs, which are CMAQ-ready emissions and converts them for use in GEOS-Chem. The tutorial uses "gridded reports", which contain annual data.  To make this operational, you instead repeat the process for each month.
# 
# Steps:
# 1. install libraries
# 2. Download data from EPA
# 3. Regrid and convert to fluxes in GEOS-Chem format
# 

# # Install Libraries
# 
# System Package | Use                   | Python Package | Use
# ---------------|-----------------------| ---------------|---------
# libgeos-dev    | Geospatial processing | basemap        | Mapping
# cdo            | Conservative regrid   | PseudoNetCDF   | IOAPI-like support
# 
# 
# * Run the next three cells

import sys

month = sys.argv[1]


# # Restart the Runtime
# 
#   * choose the "Runtime" menu, and
#   * then "Restart Runtime"
# 

# # Import Libraries and Set Earth Radius
# 
# * Importing necessary libraries
# * Setting IOAPI earth radius to WRF radius

# In[1]:


import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import PseudoNetCDF as pnc
from calendar import monthrange
import pandas as pd
os.environ['IOAPI_ISPH'] = '6370000.'


# # Download Files from EPA's newftp
# 
# * Downloading inline and gridded anthropogenic emissions
# * Downloading biogenic and non-US sectors separately
# 
# Note: if failing, remove `-q` to see errros.

# In[2]:


areaonly_sectors = sectors = [
  'afdust_adj', 'ag', 'airports',
  'nonpt', 'nonroad', 'np_oilgas', 'onroad', 'onroad_ca_adj', 'rail', 'rwc', # domestic
]
intl_sectors = [
  'onroad_can', 'onroad_mex', 'othafdust_adj', 'othar', 'othptdust_adj', 'emln_othpt', # other countries (Can/Mex)
]

natural_sectors = [
  'ocean_cl2', 'beis',
  'emln_ptagfire', 'emln_ptfire', 'emln_ptfire_othna' # fires include some anthro
]
merged_sectors = [
  'mrggrid_nobeis_norwc', # sum of areaonly_sectors and allinln - partinln sectors
]

ptonly_sectors = ['emln_cmv_c1c2_12', 'emln_cmv_c3_12', 'emln_ptegu']

allinln_sectors = [
  'emln_pt_oilgas_allinln',
  'emln_ptnonipm_allinln',
]

partinln_sectors = [
  'emln_pt_oilgas',
  'emln_ptnonipm',
]

download_sectors = areaonly_sectors + ptonly_sectors + allinln_sectors + partinln_sectors + merged_sectors + intl_sectors + natural_sectors

# Not including natural because it would duplicate
# GEOS-Chem MEGAN, ocean, and fires
# Should cmv_c1_c2 be in? cmv_c3?
# airports?
include_sectors = areaonly_sectors + ptonly_sectors + allinln_sectors + intl_sectors


# In[3]:


#month = '02'
mon3c = pd.to_datetime('2016-{}-01'.format(month)).strftime('%b').lower()
filetmpl = "{dir}/2016fh_16j_{{sector}}_{dom}_month_{{month}}.ncf".format
smoketmpl = filetmpl(dir='input', dom='12US1').format
hemcotmpl = filetmpl(dir='output', dom='0pt1degree').format
urltmpl = filetmpl(dir="ftp://newftp.epa.gov/Air/emismod/2016/v1/gridded/monthly_netCDF/", dom='12US1').format


# In[ ]:


os.makedirs('input', exist_ok=True)
for sector in download_sectors:
  outpath = smoketmpl(sector=sector, month=month)
  url = urltmpl(sector=sector, month=month)
  print(url)
  os.system(f"wget -O {outpath}  --continue -q {url}")

os.system(f"wget --continue -q ftp://newftp.epa.gov/aqmg/global/gadm/gadm36_12US1.IOAPI.nc")


# # Opening Files For Reading And Plotting
# 

# In[ ]:


smokepaths = {
    sector: smoketmpl(sector=sector, month=month)
    for sector in include_sectors + natural_sectors
}
smokefiles = {
    sector: pnc.pncopen(path, format='ioapi', mode='r')
    for sector, path in smokepaths.items()
}

reffile = smokefiles[include_sectors[0]]


# # Store Grid Parameters for Later Use
# 
# * Regridding requires knowing about the grid structure
# * We are pulling all the metadata, so that we can use what we need.
# 

# In[ ]:


gridproperties = reffile.getncatts()
exec('nominalarea = XCELL * YCELL', None, gridproperties)
exec('false_easting = -XORIG', None, gridproperties)
exec('false_northing = -YORIG', None, gridproperties)
exec('gridsize = NCOLS * NROWS', None, gridproperties)


# # Creating Custom Merge files
# 
# * Disabling warnings
#   * warns that  variables missing in the right hand file
#   * and are simply copied from the first file.
# * We don't need to see that.

# In[ ]:


def combine(f1, f2):
  """
  Arguments
  ---------
  f1 : netcdf-like object
    left-hand side of combine
  f2 : netcdf-like object
    right-hand side of combine

  Returns
  -------
  outf : netcdf-like object
    has all variables from f2 added to f1. If a variable is in f2 but
    not f1, it is added
  
  Notes
  -----
  Warnings are about variables in f1 that are not in f2, so they are ignored
  """
  import warnings
  warnings.filterwarnings('ignore')

  outf = f1 + f2
  for k, v in f2.variables.items():
      if k not in outf.variables:
          outf.copyVariable(v, key=k)
  warnings.resetwarnings()
  return outf

def merge(files):
  """
  Arguments
  ---------
  files : iterable of netcdf-like files
    files to be combined

  Returns
  -------
  outf : netcdf-like object
    has all mass from all variables in all files
  """
  outf = files[0].copy()
  for othf in files[1:]:
    outf = combine(outf, othf)
    
  return outf

def merge_from_dict(source, keys, outpath):
  """
  Arguments
  ---------
  source : dictionary
    contains netcdf-like variables keyed by sector
  keys : iterable
    keys to use to get netcdf-like variables for merging
  outpath : str
    save the merged file to disk
  
  Returns
  -------
  None
  """
  if os.path.exists(outpath):
    os.remove(outpath)
  outf = merge([source[k] for k in keys])
  history = getattr(outf, 'HISTORY', '')
  history += '; ' + ' + '.join(keys)
  setattr(outf, 'HISTORY', history)
  ondisk = outf.save(outpath, verbose=0);
  ondisk.close()


# In[ ]:


outpath = smoketmpl(sector='merge_areaonly', month=month)
merge_from_dict(smokefiles, areaonly_sectors, outpath)

outpath = smoketmpl(sector='merge_point', month=month)
merge_from_dict(smokefiles, allinln_sectors + ptonly_sectors, outpath)

outpath = smoketmpl(sector='merge_intl', month=month)
merge_from_dict(smokefiles, intl_sectors, outpath)

outpath = smoketmpl(sector='merge_natural', month=month)
merge_from_dict(smokefiles, natural_sectors, outpath)

mergedfiles = {
    sector: pnc.pncopen(smoketmpl(sector=sector, month=month), format='ioapi', mode='r')
    for sector in ['merge_areaonly', 'merge_point', 'merge_intl', 'merge_natural']
}
outpath = smoketmpl(sector='merge', month=month)
merge_from_dict(mergedfiles, ['merge_areaonly', 'merge_point', 'merge_intl'], outpath)
mergedfiles['merge'] = pnc.pncopen(outpath, format='ioapi', mode='r')


# # Create a Database of Sums
# 
# This will be used later for checking regridded data

# In[ ]:


checkspc = 'NOX'
originalrate = {
  sector: smokefiles[sector].variables[checkspc][:].sum()
  for sector in include_sectors if checkspc in smokefiles[sector].variables
}
mergedrate = {
    k: v.variables[checkspc][:].sum() for k, v in mergedfiles.items()
}


# In[ ]:


originalrate


# In[ ]:


print(sum(originalrate.values()))
print(mergedrate)


# # Comparing to US Report
# 
# * Download reference report for comparison
# * This report contains US-only sectors
#   * for cmv and beis, the smoke files include areas outside the US
#   * The smoke files should always include more mass than the report.
#   * All other sectors in the monthly report should match.
# * Comparing US sectors excluding CMV, beis, and fire for total

# In[ ]:


os.system("wget --continue ftp://newftp.epa.gov/Air/emismod/2016/v1/reports/2016fh_county_monthly_report_22jan2020_v0.csv")


# In[ ]:


monthlycounty = pd.read_csv('2016fh_county_monthly_report_22jan2020_v0.csv', comment='#')


# In[ ]:


monthlytotal = monthlycounty.groupby(['country_cd', 'census_tract_cd', 'poll']).sum().filter(like='_value').round(0)
spcmonthlytotal = monthlytotal.query('poll == "{}"'.format(checkspc)).copy()
for sector, rate in originalrate.items():
  sectork = sector.replace('emln_', '').replace('_allinln', '')
  key = ('US', sectork, checkspc)
  if key in spcmonthlytotal.index:
    spcmonthlytotal.loc[key, 'check'] = rate
  else:
    print('Skipping ', key, 'not in report')
spcmonthlytotal.round(0)


# In[ ]:


reportmatch = (
    mergedrate['merge']
    - mergedrate['merge_intl'] - originalrate['emln_cmv_c1c2_12'] - originalrate['emln_cmv_c3_12']
)
reporttotal = spcmonthlytotal.query('census_tract_cd not in ("beis", "ptfire", "ptagfire", "cmv_c1c2_12", "cmv_c3_12")').sum()
print('Report without fire, beis, cmv')
print(reporttotal.round(0))
print ('\nFor comparing to report', reportmatch)
print('Ratio', reportmatch / reporttotal[mon3c + '_value'])


# # Prepare for Regridding By Defining Grids
# 
# * EPA12US1.grid defines a Lambert Conic Conformal grid
# * latlon.grid defines a regional latitude/longitude grid
# * z.grid defines the first layer of GEOS-Chem's hybrid coordinate
# 

# In[ ]:


lccpath = 'EPA12US1.grid'
llpath = 'latlon.grid'
with open(lccpath, mode='w') as gf:
  gf.write("""gridtype = projection
gridsize = {gridsize}
xsize = {NCOLS}
ysize = {NROWS}
xinc = {XCELL}
yinc = {YCELL}
xname = x
yname = y
grid_mapping_name = lambert_conformal_conic
latitude_of_projection_origin = {YCENT}
longitude_of_projection_origin = {XCENT}
longitude_of_central_meridian = {XCENT}
standard_parallel = {P_ALP}, {P_BET}
earth_radius = 6370000.
semi_major_axis = 6370000.
semi_minor_axis = 6370000.
false_easting = {false_easting}
false_northing = {false_northing}
""".format(**gridproperties))
  
with open(llpath, mode='w') as gf:
  gf.write("""gridtype  = lonlat
gridsize  = 360000
xsize     = 900
ysize     = 400
xname     = lon
xlongname = "longitude"
xunits    = "degrees_east"
yname     = lat
ylongname = "latitude"
yunits    = "degrees_north"
xfirst    = -139.95
xinc      = 0.1
yfirst    = 20.05
yinc      = 0.1
""")
  
with open('z1.grid', mode='w') as gf:
  gf.write("""zaxistype = hybrid 
size = 1
levels = 0
positive = "up"
long_name = "GEOS-Chem levels"
axis = "Z"
""")


# # Initialize the CDO Library and Prepare Command
# 
# * Initialization is easy
# * You may get a warning  
# 

# In[ ]:


from cdo import Cdo
cdo = Cdo()
# cdo.CDO = '/work/ROMO/anaconda3/envs/geo/bin/cdo


# # Now Configure the CDO Command
# 
# cmd has 5 parts (run and described in bottom-to-top order)
# 
# * Set the vertical axis,
# * rename TSTEP axis to time, remove TFLAG, set the time
# * convert from tons/year to tons/s, tons/m2/s, then kg/m2/s
# * Set units to mass fluxes
# * Set grid to LCC and Regrid to the lon/lat grid
# 

# In[ ]:


cmdtmpl = (
    '-remapycon,{llpath} -setgrid,{lccpath} ' + # assign LCC grid and regrid
    '-setattribute,NOX@units=kgNO2/m2/s -setattribute,NO@units=kgNO2/m2/s -setattribute,*@units=kg/m2/s ' + # fix units
    '-mulc,{unitconversion} -divc,{nominalarea} -divc,{seconds} ' + # convert tons/year to kg/m2/s
    '-setreftime,1970-01-01,00:00:00,hours -settaxis,2016-{month}-01,00:00:00 -delname,TFLAG -chname,TSTEP,time ' + # fix time properties
    '-setzaxis,z1.grid -chname,LAY,lev ' +
    '{inpath}' # fix dimensions
)


# # Regrid Each file
# 
# * Using cmd described above and applying 0 to missing values
# * Output
#   * NetCDF4
#   * grid chunking
#   * level 1 compressions
# 

# In[ ]:


remakeall = False
fileopts = '-O -f nc4 -k grid -z zip_1'
os.makedirs('output', exist_ok=True)
options = dict(
  llpath=llpath, lccpath=lccpath, nominalarea=gridproperties['nominalarea'],
  unitconversion=907.185, month=month
)

cdo_sectors = include_sectors + list(mergedfiles)

for sector in cdo_sectors:
  inpath = smoketmpl(sector=sector, month=month)
  outpath = hemcotmpl(sector=sector, month=month)
  if os.path.exists(outpath) and not remakeall:
    print('Keeping cached:', outpath)
    continue
  else:
    print(outpath, end='', flush=True)
  f = pnc.pncopen(inpath, format='ioapi')
  sdate = f.getTimes()[0]
  ndays = monthrange(sdate.year, sdate.month)[1]
  options['seconds'] = ndays * 24 * 3600
  if os.path.exists(outpath):
    os.remove(outpath)
  cmd = cmdtmpl.format(inpath=inpath, **options)
  print(cmd,end='.')
  cdo.setmisstoc(0, input=cmd, output=outpath, options=fileopts)
  print('.done', flush=True)
print('Success: all complete')


# # Check Mass Conservation
# 
# * Calculate the Mass from Regridded Fluxes
# * Get the Area of the New Grid
# * Convert to mass comparable to CMAQ
#   * Multiply the area by the fluxes
#   * Convert per second rates to per year
#   * Convert kg to short tons
# 
# 

# In[ ]:


areaf = cdo.gridarea(llpath, input='-remapnn,%s -stdatm,0' % llpath, options='-f nc4', returnCdf=True)
area = areaf.variables['cell_area'][:]


# In[ ]:


hemco2dfiles = {
  sector: pnc.pncopen(hemcotmpl(sector=sector, month=month), format='netcdf', mode='r')
  for sector in include_sectors
}

hemco2dmergedfiles = {
  sector: pnc.pncopen(hemcotmpl(sector=sector, month=month), format='netcdf', mode='r')
  for sector in list(mergedfiles)
}
hemco2drate = {}
for sector in include_sectors:
  sfile = hemco2dfiles[sector]
  if checkspc not in sfile.variables:
    continue
  sdate =  sfile.getTimes()[0]
  ndays = monthrange(sdate.year, sdate.month)[1]
  checkval = (sfile.variables[checkspc][:].sum((0, 1)) * area).sum() * ndays * 24 * 3600 * 0.00110231
  hemco2drate[sector] = checkval

hemco2dmergedrate = {}
for sector in list(mergedfiles):
  sfile = hemco2dmergedfiles[sector]
  if checkspc not in sfile.variables:
    continue
  sdate =  sfile.getTimes()[0]
  ndays = monthrange(sdate.year, sdate.month)[1]
  checkval = (sfile.variables[checkspc][:].sum((0, 1)) * area).sum() * ndays * 24 * 3600 * 0.00110231
  hemco2dmergedrate[sector] = checkval


# In[ ]:


hemco2dtotal = sum(hemco2drate.values())
hemco2dmergedtotal = (hemco2dmergedrate['merge_areaonly'] + hemco2dmergedrate['merge_point'] + hemco2dmergedrate['merge_intl'])
hemco2dreportmatch = (
    hemco2dmergedrate['merge']
    - hemco2dmergedrate['merge_intl'] - hemco2drate['emln_cmv_c1c2_12'] - hemco2drate['emln_cmv_c3_12']
)
print('Report without fire, beis, cmv')
print(reporttotal.round(0))
print ('\nFor comparing to report', hemco2dreportmatch)
print('Ratio', hemco2dreportmatch / reporttotal[mon3c + '_value'])

hemco2dtotal = sum(hemco2drate.values())
hemco2dmergedtotal = (hemco2dmergedrate['merge_areaonly'] + hemco2dmergedrate['merge_point'] + hemco2dmergedrate['merge_intl'])
print('ORIG', mergedrate['merge'])
print('HEMCO2D merged ar/pt/intl', hemco2dtotal)
print('HEMCO2D single merged', hemco2dmergedtotal)
print('ORIG/HEMCO2D', hemco2dtotal / mergedrate['merge'])
hemco2dmergedrate


# In[ ]:


for key in 'NOX APIN BPIN ISOP SESQ TERP'.split():
  for sector in include_sectors:
    sfile = hemco2dfiles[sector]
    if key not in sfile.variables:
        continue
    vals = sfile.variables[key][:]
    print('{}: {} sum={:.7g}'.format(sector, key, (vals * area).sum() * 3600 * 24 * 366 * 0.00110231))


# # Plot Results
# 
# * Always good to take a look at NOX
# * Showing 2D results
# * cmv c1, c2, c3 has international emissions, which explains the small emission masses in Canada.
# 


# # Vertical Allocation Table
# 
# * Oversimplifying vertical allocation
#   * Similar to Simpson et al. 2003 Table 4.1, we apply a single representative allocation to each sector.
#   * Simpson D., Fagerli, H., Jonson, J. E., Tsyro, S., Wind, P., and Tuovinen, J.: "Transboundary Acidification, Eutrophication, and Ground Level Ozone in Europe - Part I: Unified EMEP Model Description”, EMEP Status Report 2003, The  Norwegian Meteorological Institute, Oslo, 25 Norway, 2003. https://www.emep.int/publ/reports/2003/emep_report_1_part1_2003.pdf

# In[ ]:


import pandas as pd
import io

vertical_allocation = pd.read_csv(io.StringIO("""L,bottom_m,top_m,othpt,ptegu,ptnonipm,pt_oilgas,ptfire,cmv_c3
1,-6,123,0.806067059,0.027331187,0.751204553,0.918700528,0.11976381,0.537538889
2,123,254,0.087955253,0.12260809,0.144224191,0.059678397,0.098297847,0.338982929
3,254,387,0.053252808,0.192017292,0.052822078,0.013323237,0.08277491,0.123478182
4,387,521,0.025605531,0.212736741,0.023840537,0.003978944,0.074659995,0
5,521,657,0.011408677,0.167549199,0.011000364,0.001724366,0.077561337,0
6,657,795,0.005787852,0.110491495,0.005661267,0.001038275,0.073349589,0
7,795,934,0.003758579,0.063653198,0.003165207,0.000633395,0.066694526,0
8,934,1075,0.002673887,0.036979211,0.00197344,0.000387233,0.061845442,0
9,1075,1218,0.001675656,0.022844579,0.001300664,0.000243355,0.058324971,0
10,1218,1363,0.000670159,0.013909341,0.000861865,0.000122994,0.053902287,0
11,1363,1510,0.000357644,0.008781226,0.000621856,6.48583E-05,0.044774262,0
12,1510,1659,0.00022693,0.006647819,0.00052355,4.04985E-05,0.041122289,0
13,1659,1860,0.000159139,0.004305749,0.000466642,2.12041E-05,0.0306914,0
14,1860,2118,0.000152154,0.003380562,0.000680222,1.75333E-05,0.031529605,0
15,2118,2382,0.000113674,0.001827346,0.000612964,1.03278E-05,0.02551535,0
16,2382,2654,8.79442E-05,0.000957313,0.000260885,5.76406E-06,0.020748307,0
17,2654,2932,3.14778E-05,0.000523901,0.000311535,3.23215E-06,0.014366484,0
18,2932,3219,1.29957E-05,0.000381747,0.000231719,2.27428E-06,0.010906613,0
19,3219,3665,1.37334E-06,0.000409744,0.000115745,1.78044E-06,0.009224629,0
20,3665,4132,3.15633E-07,0.000340645,4.19637E-05,6.00909E-07,0.003409241,0
21,4132,4623,3.71912E-07,0.00026689,1.84353E-05,4.14262E-07,0.00014894,0
22,4623,5142,4.01952E-07,0.000231509,1.03385E-05,3.73296E-07,0.000125351,0
23,5142,5692,1.17457E-07,0.000225471,9.11318E-06,3.01765E-07,8.76892E-05,0
24,5692,6277,0,0.000310076,9.82452E-06,1.13299E-07,3.24343E-05,0
25,6277,6905,0,0.000327088,1.07925E-05,0,3.0049E-05,0
26,6905,7582,0,0.000260516,1.2311E-05,0,2.85539E-05,0
27,7582,8320,0,0.000130502,4.58591E-06,0,2.08973E-05,0
28,8320,9409,0,0.000119402,2.72213E-06,0,3.26732E-05,0
29,9409,10504,0,0.000106184,6.30588E-07,0,1.74394E-05,0
30,10504,11578,0,9.07766E-05,0,0,1.30783E-05,0
31,11578,12633,0,0.000117164,0,0,0,0
32,12633,13674,0,0.000122616,0,0,0,0
33,13674,14706,0,1.54213E-05,0,0,0,0
34,14706,15731,0,0,0,0,0,0
35,15731,16753,0,0,0,0,0,0
36,16753,17773,0,0,0,0,0,0
37,17773,19855,0,0,0,0,0,0
38,19855,22004,0,0,0,0,0,0
39,22004,24240,0,0,0,0,0,0
40,24240,26596,0,0,0,0,0,0
41,26596,31716,0,0,0,0,0,0
42,31716,37574,0,0,0,0,0,0
43,37574,44286,0,0,0,0,0,0
44,44286,51788,0,0,0,0,0,0
45,51788,59924,0,0,0,0,0,0
46,59924,68392,0,0,0,0,0,0
47,68392,80581,0,0,0,0,0,0
""")).set_index(['L'])


# # Set a Maximum Number of Layers
# 
# * Some layers have non-zero values up to layer 33
#   * This seems unreasonable.
#   * Simpson et al. 2003 had up to layer 11 in this vertical coordinate
# * Setting a maximum top layer of 11
#   * renormalizing mass within those 11 layers.

# In[ ]:


nl = 11
vertical_allocation_layers = vertical_allocation.loc[1:nl]
vertical_allocation_layer_fraction = vertical_allocation_layers / vertical_allocation_layers.sum()


# # Apply Vertical Allocation
# 
# * Map each smoke sector to a vertical allocation
# * Manually fixing meta data
#   * (hyai, hybi, and lev)
#   * renaming dimensions to ilev

# In[ ]:


import gc


sector2verticalallocation = {
  'emln_cmv_c1c2_12': 'cmv_c3',
  'emln_cmv_c3_12': 'cmv_c3',
  'emln_othpt': 'othpt',
  'emln_pt_oilgas_allinln': 'pt_oilgas',
  'emln_ptegu': 'ptegu',
  'emln_ptnonipm_allinln': 'ptnonipm'
}

for sector in ptonly_sectors + allinln_sectors:
  gc.collect()
  inpath = hemcotmpl(sector=sector, month=month)
  outpath = inpath.replace('0pt1degree', '0pt1degree_3D')
  origlayerfractions = vertical_allocation_layer_fraction[sector2verticalallocation[sector]].values
  nl = (origlayerfractions == 0).argmax()
  if nl == 0:
    nl = len(origlayerfractions)
  layerfractions = origlayerfractions[:nl]
  layerfractions /= layerfractions.sum()
  print(sector, origlayerfractions)
  print(sector, nl, layerfractions)
  print('Start.', outpath, end='')
  if os.path.exists(outpath):
    os.remove(outpath)
  emln2df = pnc.pncopen(inpath, format='netcdf')
  emln3df = emln2df.sliceDimensions(lev=[0] * nl).sliceDimensions(lev=[0] * nl)#.sliceDimensions(nhyi=[0] * (nl + 1))
  for key, var in emln3df.variables.items():
    if 'lev' in var.dimensions:
      newshape = tuple([(slice(None) if dk == 'lev' else None) for dk in var.dimensions])
      var[:] = var[:] * layerfractions[newshape]
  # hyai = emln3df.variables['hyai'][:] = [
  #   0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
  #   3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01
  # ][:nl+1]
  # hybi = emln3df.variables['hybi'][:] = [
  #   1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
  #   8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01
  # ][:nl+1]
  # emln3df.variables['hybm'][:] = np.convolve([.5, .5], hybi, mode='valid')
  # emln3df.variables['hyam'][:] = np.convolve([.5, .5], hyai, mode='valid')
  lev = emln3df.createVariable('lev', 'd', ('lev',))
  lev.long_name = "GEOS-Chem levels" ;
  lev.standard_name = "GEOS-Chem levels" ;
  lev.units = "model level" ;
  lev.positive = "up" ;
  lev[:] = np.arange(nl, dtype='i')
  # Requires special treatment in variable definition
  # to match geos-chem optimized chunking.
  outnc = emln3df.copy(variables=False).save(outpath, verbose=0)
  dimlens = {k: len(v) for k, v in emln3df.dimensions.items()}
  chunksizes = dimlens.copy()
  chunksizes['lev'] = 1
  chunksizes['time'] = 1
  for vark in emln3df.dimensions:
    if vark not in emln3df.variables:
        continue
    var = emln3df.variables[vark]
    opts = dict(complevel=0, chunksizes=tuple([dimlens[dk] for dk in var.dimensions]))
    ovar = outnc.createVariable(vark, var.dtype.char, var.dimensions, **opts)
    ovar.setncatts(var.getncatts())
    ovar[:] = var[:]
    
  for vark, var in emln3df.variables.items():
    if vark in emln3df.dimensions:
        continue
    opts = dict(complevel=1, chunksizes=tuple([chunksizes[dk] for dk in var.dimensions]))
    ovar = outnc.createVariable(vark, var.dtype.char, var.dimensions, **opts)
    ovar.setncatts(var.getncatts())
    ovar[:] = var[:]
    
  outnc.close()
  print('.Complete')


# # Now Do a Mass Conservation Check
# 
# * Open the 3d file
# * Open the reference file
# * Sum the values and compare

# In[ ]:


from glob import glob

emln3dpaths = sorted(glob('output/*_3D_*.ncf'))
for inpath in emln3dpaths:
  refpath = inpath.replace('_3D', '')
  f3d = pnc.pncopen(inpath)
  fref = pnc.pncopen(refpath)
  print(refpath)
  print('emln   {:2d} lev NOX {:.7e} sum kgNO2/m2/s'.format(len(fref.dimensions['lev']), fref.variables['NOX'][:].sum()))
  print('emln3d {:2d} lev NOX {:.7e} sum kgNO2/m2/s'.format(len(f3d.dimensions['lev']), f3d.variables['NOX'][:].sum()))
  print('Within numerical precision', np.allclose(fref.variables['NOX'][:, 0], f3d.variables['NOX'][:].sum(1)))


# # Configure Translation from CMAQ to GEOS-Chem
# 
# * Translation is not always 1:1
# * Apply scaling factors
#   * Diurnal scale factors (25, 26)
#   * Day of week (210-222)
#   * Split OC/EC (70-73)
#   * Convert NO kgNO2/m2/s to kgNO/m2/s (115)

# In[ ]:


cq2gc = {
    'ACET': [['ACET', '26/213/1007']],
    'ALD2': [['ALD2', '26/213/1007']],
    'ALDX': [['RCHO', '26/213/1007']],
    'ACROLEIN': [['MACR', '26/213/1007']],
    'BENZ': [['BENZ', '26/213/1007']],
    'ECH4': [['CH4', '1007']],
    'CL': [['Cl', '1007']],
    'CL2': [['Cl2', '1007']],
    'CLO': [['ClO', '1007']],
    'CLNO2': [['ClNO2', '1007']],
    'CLNO3': [['ClNO3', '1007']],
    'CO': [['CO', '26/211/1007']],
    'ETOH': [['EOH', '26/213/1007']],
    'ETHA': [['C2H6', '26/217/1007']],
    'ETH': [['C2H4', '26/213/1007']],
    'FACD': [['HCOOH', '26/213/1007']],
    'FORM': [['CH2O', '26/213/1007']],
    'GLY': [['GLYX', '26/214/1007']], # Following MEK
    'GLYD': [['GLYC', '26/214/1007']], # Following MEK
    'HCL': [['HCl', '1007']],
    'HOCL': [['HOCl', '1007']],
    'HONO': [['HNO2', '25/210/1007']],
    'HPLD': [['HPALD', '26/214/1007']], # Following MEK
    'IOLE': [['PRPE', '26/215/1007']],
    'KET': [['MEK', '26/214/1007']],
    'MEPX': [['MP', '26/213/1007']], # Following MEOH
    'MEOH': [['MOH', '26/213/1007']],
    'NH3': [['NH3', '26/213/1007']],
    'NO': [['NO', '115/25/210/1007']],
    'NO2': [['NO2', '25/210/1007']],
    'OLE': [['PRPE', '26/215/1007']],
    'PACD': [['MAP', '26/213/1007']], # Following MEOH
    'PAR': [['ALK4', '26/212/1007']],
    'PNA': [['HNO4', '26/213/1007']], # Following MEOH
    'PRPA': [['C3H8', '26/216/1007']],
    'SO2': [['SO2', '26/218/1007']],
    'SULF': [['SO4', '26/218/1007']],
    'TOL': [['TOLU', '26/213/1007']],
    'PEC': [['BCPI', '26/221/1007/70'], ['BCPO', '26/221/256/1007/71']],
    'PFE': [['pFe', '26/219/1007']],
    'POC': [['OCPI', '26/222/1007/72'], ['OCPO', '26/222/257/1007/73']],
    'PNH4': [['NH4', '26/218/1007']], # need to confirm this
    'PNO3': [['NIT', '26/218/1007']], # need to confirm this
    'PSO4': [['SO4', '26/219/1007']],
    'XYLMN': [['XYLE', '26/213/1007']],
}
# Ignore special species
exec('ALD2_PRIMARY = FORM_PRIMARY = NH3_FERT = []', None, cq2gc)
# ignore inventory meta variables
exec('HFLUX = VOC_INV = VOC_BEIS = CO2_INV = N2O_INV = CH4_INV = []', None, cq2gc)
# Ignore duplicative species
# TOLU in TOL
# NOX = NO + NO2 + HONO
# CH4 in ECH4
exec('CH4 = TOLU = NOX = []', None, cq2gc)
# Ignore coordinate variables
exec('hybm = hyam = []', None, cq2gc)
# Ignore APIN, BPIN, and SESQ - they were from BEIS
# and were removed
exec('APIN = BPIN = SESQ = []', None, cq2gc)
# ISOP and TERP are mostly biogenic, so igoring remainder
exec('ISOP = TERP = []', None, cq2gc)
# Ignore species GC doesn't have
exec('BUTADIENE13 = ETHY = NAPH = []', None, cq2gc)
# unreactive and nonvolatile are not used
exec('UNK = UNR = NVOL = NR = []', None, cq2gc)
# Ignore unused particulate species
exec('SOAALK = PAL = PCA = PCL = PFE = PH2O = PK = PSI = PTI = []', None, cq2gc)
exec('PM2_5 = PMC = PMG = PMN = PMOTHR = PNCOM = P = PCA = PCL = PFE = PH2O = PK = []', None, cq2gc)


# # Define HEMCO_Config.rc Writer
# 
# * Create a HEMCO_Config.rc that reads all included files
# * Only needed to be made once and editted to use $MM instead of 01..12
# 

# In[ ]:


def changepathtopattern(path):
  return '$ROOT/NEI2016fh/' + path.replace('month_' + month, 'month_$MM')

def getsector(hcpath):
    sector = hcpath.split('16j_')[1].split('_0pt')[0]
    sector = sector.replace('emln_', '').replace('_allinln', '').replace('_12', '').replace('_adj', '').replace('cmv_', '')
    sector = sector.replace('_oilgas', 'og')
    return sector

def writeconfig(sectors, outpath):
    hcpaths = []
    for sector in sectors:
      hemco2dpath = hemcotmpl(sector=sector, month=month)
      hemco3dpath = hemco2dpath.replace('0pt1degree', '0pt1degree_3D')
      if os.path.exists(hemco3dpath):
          hcpaths.append(hemco3dpath)
      elif os.path.exists(hemco2dpath):
          hcpaths.append(hemco2dpath)
      else:
          raise KeyError('Could not find regridded: ' + hemco2dpath)

    defaults = set()
    ignores = set()
    with open(outpath, 'w') as hcf:
        hcf.write('(((EPA2016_MONMEAN\n')
        for hcpath in hcpaths:
            hcpatt = changepathtopattern(hcpath)
            sector = getsector(hcpath)
            print(sector, hcpatt, end='', flush=True)
            hcfile = pnc.pncopen(hcpath, format='netcdf')
            for cqkey, v in hcfile.variables.items():
                if cqkey in hcfile.dimensions or cqkey in ('hyai', 'hybi'):
                  continue
                elif cqkey in ('TOLU',):
                  warn('TOLU mass is duplicated by TOL')
                if cqkey in cq2gc:
                  gctrans = cq2gc.get(cqkey)
                  if len(gctrans) == 0:
                    ignores.add(cqkey)
                else:
                  defaults.add(cqkey)
                  gctrans = [[cqkey, '1007']]
                for gckey, scale in gctrans:
                  if gckey in ['ACET', 'MEK', 'ALD2', 'PRPE', 'PRPA', 'BENZ', 'TOLU', 'XYLE', 'EOH', 'ALK4', 'ISOP']:
                    units = 'kgC/m2/s'
                  else:
                    units = v.units.strip()
                  opts = dict(
                    unit=units,
                    gckey=gckey,
                    cqkey=cqkey,
                    sector=sector,
                    path=hcpatt,
                    scale=scale,
                    cat='1/2',
                    hier=50
                  )
                  hcf.write('0 EPA16_{gckey}__{sector}{cqkey} {path}  {cqkey}       2016-2016/1-12/1/0 C xyz  {unit}  {gckey}   {scale}     {cat} {hier}\n'.format(**opts))
                  # If I use - to repeat the file, the mass is from the previous cqkey too.
                  # hcpatt = '-'
            print()
        hcf.write(')))EPA2016_MONMEAN\n')
    print('Ignored', sorted(ignores))
    print('Defaults', sorted(defaults))


# # Choose paths to use in HEMCO
# 
# * hemco3d files take precedence
# * hemco2d are the default
# * Must have one or the other

# In[ ]:


writeconfig(include_sectors, 'output/EPA_Config.rc')
writeconfig(['merge_areaonly', 'merge_intl'] + ptonly_sectors + allinln_sectors, 'output/EPAMERGED_Config.rc')

print('Finished')
# # Review the HEMCO file



# # Summary
# 
# You've created 3 important in the output folder.
# 1. A 3D emln file with no fires.
# 2. A 2D gridded file with no beis
# 3. A HEMCO_Config.rc file
# 
# The are two additional files that are just references:
# 1. tot_nobeis_nofire is a refernce fo the expected total
# 2. emln_nofire (non-3D)
# 
# Download the output directory

# In[ ]:


#os.system('zip output.zip output/*')



# In[ ]:




