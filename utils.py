import os
import datetime
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import CIRS, GCRS
from astropy.coordinates import SkyCoord, Distance
from astropy.table import Table
import pickle
from astroquery.xmatch import XMatch
#from astroquery.gaia import Gaia
#Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source"
#Gaia.ROW_LIMIT = 2
 
__all__ = ['modulepath', 'load_obd_module', 'write_obd', 'coord_offset', 'cal_offset', 
           'cal_coord_motion', 'sc_offset', 'get_coord_plain', 'get_coord_colon',
           'get_pos_current', 'get_pos_J2000', 'get_pos_ao', 'read_coordinate', 
           'coordinate_convert_epoch', 'search_gaia_single', 'search_gaia_2table', 
           'xmatch_gaia', 'color_GrpH', 'color_GrpK']

#-> Obtain the current path
pathList = os.path.abspath(__file__).split("/")
#-> Create the path to the gData module
modulepath = "/".join(pathList[0:-1])


def load_obd_module(tname, tpath=modulepath):
    '''
    Load the template parameters.
    
    Parameters
    ----------
    tname : string
        The obd module name.
    tpath : string (default: module path)
        The path to the obd modules.
    '''
    f = open('{0}/tpl_modules/{1}.obd'.format(tpath, tname), 'r')
    text = f.readlines()

    cDict = {}
    for loop, t in enumerate(text):

        # Skip the comments
        if t[0] == '#':
            continue

        c0 = t.split(';')

        # Skip the empty lines
        if len(c0[0]) == 1:
            continue

        c1 = c0[0].split()
        key = c1[0]
        if len(c1) == 1:
            content = ''
        else:
            content = ' '.join(c1[1:])
            c2 = content.split('#')
            content = c2[0].strip()  # Remove the white space

            if content[0] == '"':
                content = content[1:]
            if content[-1] == '"':
                content = content[:-1]

        cDict[key] = content
    return cDict


def write_obd(template_name, obd_dict, nchars=50):
    '''
    Write down the OBD.
    
    Parameters
    ----------
    template_name : string
        The full path of the OBD.
    obd_dict : dict
        The dict of all OBD modules.
    nchars : int
        The total number of characters in one line.
    '''
    # Put everything in strings
    lines = []
    for mn in obd_dict:
        lines.append('# {}\n'.format(mn))
        for key in obd_dict[mn]:
            if key == 'PAF.HDR.START':
                content = ''
            else:
                content = '"{}"'.format(obd_dict[mn][key])
            nc_left = nchars-len(key)
            if nc_left < 0:
                raise KeyError('The key ({0}) is longer than {1}!'.format(key, nchars))
            lineformat = '{{0}} {{1:>{0}}};\n'.format(nc_left)
            lines.append(lineformat.format(key, content))
        lines.append('\n')

    f = open(template_name, 'w')
    f.writelines(lines)
    f.close()


def cal_coord_motion(c, pma=None, pmd=None, plx=None, radvel=None,
                     time_ref='J2000.0', time_cal='now', frame=GCRS):
    """
    Calculate the current coordinates considering the motion of the source.

    Parameters
    ----------
    c : Astropy SkyCoord
        The coordinate to calculate the motion.
    pma (optional) : Astropy Quantity (angular velocity)
        The proper motion of RA.
    pmd (optional) : Astropy Quantity (angular velocity)
        The proper motion of DEC.
    plx (optional) : Astropy Quantity (angle)
        The parallex.
    radvel (optional) : Astropy Quantity (velocity)
        The radial velocity.
    time_ref : string (default: 'J2000.0')
        The reference time of the input data.
    time_cal : string (default: 'now')
        The time to calculate the motion of the coordinate.
    frame : coordinate systems (default: GCRS)
        The coordinate system, which need to consider the orbit of the Earth
        around the Sun.

    Returns
    -------
    c_c : Astropy SkyCoord
        The calculated coordinate.
    """
    if pma is None:
        pma = 1e-16*u.mas/u.yr
    if pmd is None:
        pmd = 1e-16*u.mas/u.yr
    if plx is None:
        plx = 0 * u.mas
    if radvel is None:
        radvel = 0 * u.km/u.s
    c_m = SkyCoord(ra=c.ra, dec=c.dec, pm_ra_cosdec=pma, pm_dec=pmd,
                   distance=Distance(parallax=plx), radial_velocity=radvel,
                   frame=c.frame.name, obstime=time_ref)
    if time_cal == 'now':
        time_cal = Time.now()
    else:
        time_cal = Time(time_cal)
    c_c = c_m.apply_space_motion(time_cal).transform_to(GCRS(obstime=time_cal))
    return c_c


def sc_offset(c_sc, c_ft, pma_sc=None, pmd_sc=None, plx_sc=None, radvel_sc=None,
              pma_ft=None, pmd_ft=None, plx_ft=None, radvel_ft=None,
              time_ref="J2000.0", time_cal='now', frame=GCRS):
    '''
    Calculate the offset of the SC fiber from the FT target.

    Parameters
    ----------
    c_sc : Astropy SkyCoord
        The coordinate of the SC target.
    c_ft : Astropy SkyCoord
        The coordinate of the FT target.
    pma_sc (optional) : float
        The proper motion of RA, units: mas/yr.
    pmd_sc (optional) : float
        The proper motion of DEC, units: mas/yr.
    plx_sc (optional) : float
        The parallex, units: mas.
    radvel_sc (optional) : float
        The radial velocity, units: km/s.
    pma_ft (optional) : float
        The proper motion of RA, units: mas/yr.
    pmd_ft (optional) : float
        The proper motion of DEC, units: mas/yr.
    plx_ft (optional) : float
        The parallex, units: mas.
    radvel_ft (optional) : float
        The radial velocity, units: km/s.

    Returns
    -------
    delta_ra : float
        The offset of RA, units: mas.
    delta_dec : float
        The offset of DEC, units: mas.
    '''
    if not pma_sc is None:
        pma_sc = pma_sc * u.mas / u.yr
    if not pmd_sc is None:
        pmd_sc = pmd_sc * u.mas / u.yr
    if not plx_sc is None:
        plx_sc = plx_sc * u.mas
    if not radvel_sc is None:
        radvel_sc = radvel_sc * u.km/u.s
    if not pma_ft is None:
        pma_ft = pma_ft * u.mas / u.yr
    if not pmd_ft is None:
        pmd_ft = pmd_ft * u.mas / u.yr
    if not plx_ft is None:
        plx_ft = plx_ft * u.mas
    if not radvel_ft is None:
        radvel_ft = radvel_ft * u.km/u.s
    c_sc_now = cal_coord_motion(c_sc, pma=pma_sc, pmd=pmd_sc, plx=plx_sc, radvel=radvel_sc,
                                time_ref=time_ref, time_cal=time_cal, frame=frame)
    c_ft_now = cal_coord_motion(c_ft, pma=pma_ft, pmd=pmd_ft, plx=plx_ft, radvel=radvel_ft,
                                time_ref=time_ref, time_cal=time_cal, frame=frame)
    delta_ra, delta_dec = cal_offset(c_sc_now, c_ft_now)
    return delta_ra, delta_dec


def coord_offset(delta_ra, delta_dec, c0, frame='icrs'):
    """
    Calculate the coordinate, that is offset from a reference coordinate.

    Parameters
    ----------
    delta_ra : float
        The ra offset of the target, units: arcsec.
    delta_dec : float
        The dec offset of the target, units: arcsec.
    c0 : Astropy SkyCoord
        The reference coordinate.
    frame : string (default: 'icrs')
        The coordinate frame.

    Returns
    -------
    c : Astropy SkyCoord
        The coordinates of the target.
    """
    ra_t  = c0.ra.arcsec + delta_ra / np.cos(c0.dec.radian)
    dec_t = c0.dec.arcsec + delta_dec
    c_t = SkyCoord(ra=ra_t*u.arcsec, dec=dec_t*u.arcsec, frame=frame)
    return c_t


def cal_offset(c, c_ref, units='mas'):
    """
    Calculate the coordinate offset betweenn the two coordinates.

    Parameters
    ----------
    c : Astropy SkyCoord
        The coordinate to calculate the offset.
    c_ref : Astropy SkyCoord
        The reference coordinate.

    Returns
    -------
    (delta_ra, delta_dec) : (float, float)
        The offsets of right ascension and declination, units: arcsec.
    """
    sep = c_ref.separation(c)
    pa = c_ref.position_angle(c)
    delta_ra = (sep * np.sin(pa)).to(units)
    delta_dec = (sep * np.cos(pa)).to(units)
    return (delta_ra, delta_dec)


def get_coord_plain(c):
    '''
    Convert the coordinate to the (HHMMSS.SSS, DDMMSS.SSS) format.
    
    Parameters
    ----------
    c : SkyCoord
        The coordinate object.
    
    Returns
    -------
    ra_plain, dec_plain : string
        The coordinates in HHMMSS.SSS and DDMMSS.SSS
    '''
    ra, dec = c.to_string('hmsdms').split(' ')
    ra_h = ra[:2]
    ra_m = ra[3:5]
    ra_s = float(ra[6:-1])
    dec_d, dec_tmp = dec.split('d')
    dec_m, dec_s = dec_tmp.split('m')
    dec_s = float(dec_s[:-1])
    ra_plain = '{0}{1}{2}'.format(ra_h, ra_m, '{0:.3f}'.format(ra_s).zfill(6))
    dec_plain = '{0}{1}{2}'.format(dec_d, dec_m, '{0:.3f}'.format(dec_s).zfill(6))
    if dec_plain[0] == '+':
        dec_plain = dec_plain[1:]
    return ra_plain, dec_plain


def get_coord_colon(c):
    '''
    Convert the coordinate to the (HH:MM:SS.SSS, DD:MM:SS.SSS) format.
    
    Parameters
    ----------
    c : SkyCoord
        The coordinate object.
    
    Returns
    -------
    ra_colon, dec_colon : string
        The coordinates in HH:MM:SS.SSS and DD:MM:SS.SSS
    '''
    try:
        ra, dec = c.to_string('hmsdms').split(' ')
    except ValueError:
        return '00:00:00.000', '00:00:00.000'
    ra_h = ra[:2]
    ra_m = ra[3:5]
    ra_s = float(ra[6:-1])
    dec_d, dec_tmp = dec.split('d')
    dec_m, dec_s = dec_tmp.split('m')
    dec_s = float(dec_s[:-1])
    ra_colon = '{0}:{1}:{2}'.format(ra_h, ra_m, '{0:.3f}'.format(ra_s).zfill(6))
    dec_colon = '{0}:{1}:{2}'.format(dec_d, dec_m, '{0:.3f}'.format(dec_s).zfill(6))
    return ra_colon, dec_colon


def get_pos_current(ra, dec, pma, pmd, parallax, radvel):
    '''
    Get the target position dict for the sequencer.
    
    Parameters
    ----------
    ra : float
        RA in degree.
    dec : float
        DEC in degree.
    pma : float
        Proper motion in mas.
    pmd : float
        Proper motion in mas.
    plx : float
        Parallax in mas.
    radvel : float
        Radial velocity in km/s.
        
    Returns
    -------
    pos : dict
        The dict of the position information of the target.
    '''
    c = read_coordinate(ra, dec)
    c_cur = cal_coord_motion(c, pma*u.mas/u.yr, pmd*u.mas/u.yr, parallax*u.mas, radvel*u.km/u.s)
    ra_hms, dec_dms = get_coord_plain(c_cur)
    pos = dict(ra=ra_hms, 
               dec=dec_dms, 
               pma='{0:.5f}'.format(pma*1e-3), 
               pmd='{0:.5f}'.format(pmd*1e-3), 
               parallax='{0:.5f}'.format(parallax*1e-3), 
               radvel=radvel)
    return pos


def get_pos_J2000(ra, dec, pma, pmd, parallax, radvel):
    '''
    Get the target position dict for the sequencer.
    
    Parameters
    ----------
    ra : float
        RA in degree.
    dec : float
        DEC in degree.
    pma : float
        Proper motion in mas.
    pmd : float
        Proper motion in mas.
    plx : float
        Parallax in mas.
    radvel : float
        Radial velocity in km/s.
        
    Returns
    -------
    pos : dict
        The dict of the position information of the target.
    '''
    c = read_coordinate(ra, dec)
    ra_hms, dec_dms = get_coord_plain(c)
    pos = dict(ra=ra_hms, 
               dec=dec_dms, 
               pma='{0:.5f}'.format(pma*1e-3), 
               pmd='{0:.5f}'.format(pmd*1e-3), 
               parallax='{0:.5f}'.format(parallax*1e-3), 
               radvel=radvel)
    return pos


def get_pos_ao(ra, dec, pma, pmd, parallax, radvel):
    '''
    Get the AO source position dict for the sequencer.
    
    Parameters
    ----------
    ra : float
        RA in degree.
    dec : float
        DEC in degree.
    pma : float
        Proper motion in mas.
    pmd : float
        Proper motion in mas.
    plx : float
        Parallax in mas.
    radvel : float
        Radial velocity in km/s.
        
    Returns
    -------
    pos : dict
        The dict of the position information of the target.
    '''
    c = read_coordinate(ra, dec)
    #c_cur = cal_coord_motion(c, pma*u.mas/u.yr, pmd*u.mas/u.yr, parallax*u.mas, radvel*u.km/u.s)
    ra_hms, dec_dms = get_coord_colon(c)
    pos = dict(ra=ra_hms, 
               dec=dec_dms, 
               pma='{0:.5f}'.format(pma*1e-3), 
               pmd='{0:.5f}'.format(pmd*1e-3), 
               parallax='{0:.5f}'.format(parallax*1e-3), 
               radvel=radvel)
    return pos
    
    
def coordinate_convert_epoch(ra, dec, pma, pmd, parallax, radial_velocity, 
                             ref_epoch, target_epoch):
    '''
    Convert the coordinate to the targeted epoch.
    
    Parameters
    ----------
    ra : float
        R.A. in degree.
    dec : float
        Decl. in degree.
    pma : float
        Proper motion of R.A., units: mas/yr.
    pmd : float
        Proper motion of Decl., units: mas/yr.
    parallax : float
        Parallax, in mas.
    radial_velocity : float
        Radial velocity, units: km/s.
    ref_epoch : string
        The reference epoch. J2015.5 for Gaia.
    target_epoch : string
        The targeted epoch. Typically, J2000.
    
    Returns
    -------
    c_t : SkyCoord
        The converted coordinate.
    '''
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, 
                 pm_ra_cosdec=pma*u.mas/u.yr, pm_dec=pmd*u.mas/u.yr, 
                 distance=Distance(parallax=parallax*u.mas), 
                 radial_velocity=radial_velocity*u.km/u.s, 
                 obstime=ref_epoch)
    c_t = c.apply_space_motion(Time(target_epoch))
    return c_t
    
    
def search_gaia_single(ra, dec, radius):
    '''
    Search the target in Gaia.
    
    Parameters
    ----------
    ra : string or float
        R.A. in 'h:m:s' or degree.
    dec : string or float
        Decl. in 'd:m:s' or degree.
    radius : float
        Searching radius in arcsec.
        
    Returns
    -------
    Result table of the closest target.
    '''
    coord = read_coordinate(ra, dec)
    radius = u.Quantity(radius, u.arcsec)
    j = Gaia.cone_search_async(coordinate=coord, radius=radius)
    r = j.get_results()
    
    len_tb = len(r)
    if len_tb == 0:
        print(ra, dec)
        raise RuntimeError('Cannot find any results!')
    
    if len_tb > 1:
        fltr = np.arange(len_tb) == np.argmin(r['dist'])
        r = r[fltr]
    
    return r
    
    
def read_coordinate(ra, dec):
    '''
    Read in the coordinate, either in degree or hourangle. Only use ICRS frame.
    
    Parameters
    ----------
    ra : float or string
        The right ascension (degree or HH:MM:SS).
    dec : float or string
        The declination (degree or DD:MM:SS).
    
    Returns
    -------
    c : SkyCoord
        The coordinate object.
    '''
    if isinstance(ra, str):
        assert isinstance(dec, str)
        c = SkyCoord('{0} {1}'.format(ra, dec), frame='icrs', unit=(u.hourangle, u.deg))
    else:
        c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    return c
    
    
def search_gaia_2table(ra, dec, radius):
    '''
    Search for Gaia info and convert them to J2000 into a table.
    '''
    r = search_gaia_single(ra, dec, radius)
    
    ra = r['ra'][0]
    dec = r['dec'][0]
    pma = r['pmra'][0]
    pmd = r['pmdec'][0]
    plx = r['parallax'][0]
    rv = r['radial_velocity'][0]
    ref_epoch = 'J{}'.format(r['ref_epoch'][0])
    target_epoch = 'J2000'
    
    if plx < 0:
        plx = 0
    c2000 = coordinate_convert_epoch(ra, dec, pma, pmd, plx, rv, ref_epoch, target_epoch)
    
    ra_c, dec_c = get_coord_colon(c2000)
    
    tb = Table([[ra_c], [dec_c], [pma], [pmd], [plx], [rv], r['phot_g_mean_mag']], 
               names=['ra', 'dec', 'pma', 'pmd', 'plx', 'rv', 'G'])
    tb['pma'].format = '%.3f'
    tb['pmd'].format = '%.3f'
    tb['plx'].format = '%.3f'
    tb['rv'].format = '%.2f'
    tb['G'].format = '%.1f'
    return tb


def xmatch_gaia(t, radius, colRA, colDec, vizier_code='I/345/gaia2'):
    '''
    Cross match the table with Gaia.
    
    Parameters
    ----------
    t : Astropy Table
        The table of targets.
    radius : float
        The cross match radius, units: arcsec.
    colRA : string
        The column name of the RA.
    colDec : string
        The column name of the Dec.
    vizier_code : string
        The vizieR code of the Gaia table.
        
    Returns
    -------
    t_f : Astropy Table
        The table of cross matched results.
    '''
    for cn in ['ra_gaia_J2000', 'dec_gaia_J2000', 'pma', 'pmd', 'plx', 'rv', 'G']:
        assert cn not in t.colnames, 'The input table has {} as a column!'.format(cn)
    
    t_o = XMatch.query(cat1=t,
                       cat2='vizier:{}'.format(vizier_code),
                       max_distance=radius * u.arcsec, colRA1=colRA,
                       colDec1=colDec, colRA2='RAJ2000', colDec2='DECJ2000')
    
    ra_j2000 = []
    dec_j2000 = []
    ra_j2015 = []
    dec_j2015 = []
    pma = []
    pmd = []
    plx = []
    rv = []
    G = []
    Grp = []
    Gbp = []
    for loop in range(len(t)):
        fltr = (t_o[colRA] == t[colRA][loop]) & (t_o[colDec] == t[colDec][loop])
        if np.sum(fltr) == 0:
            ra_j2000.append(np.nan)
            dec_j2000.append(np.nan)
            ra_j2015.append(np.nan)
            dec_j2015.append(np.nan)
            pma.append(np.nan)
            pmd.append(np.nan)
            plx.append(np.nan)
            rv.append(np.nan)
            G.append(np.nan)
            Grp.append(np.nan)
            Gbp.append(np.nan)
        else:
            t_f = t_o[fltr]
            idx = np.argmin(t_f['angDist'])
            ra_j2000.append(t_f['ra_epoch2000'][idx])
            dec_j2000.append(t_f['dec_epoch2000'][idx])
            ra_j2015.append(t_f['ra'][idx])
            dec_j2015.append(t_f['dec'][idx])
            pma.append(t_f['pmra'][idx])
            pmd.append(t_f['pmdec'][idx])
            plx.append(t_f['parallax'][idx])
            rv.append(t_f['radial_velocity'][idx])
            G.append(t_f['phot_g_mean_mag'][idx])
            Grp.append(t_f['phot_rp_mean_mag'][idx])
            Gbp.append(t_f['phot_bp_mean_mag'][idx])
    t_f = t.copy()
    t_f.add_columns([ra_j2000, dec_j2000, ra_j2015, dec_j2015, pma, pmd, plx, rv, G, Grp, Gbp], 
                   names=['ra_gaia_J2000', 'dec_gaia_J2000', 'ra_gaia_J2015', 'dec_gaia_J2015',
                          'pma', 'pmd', 'plx', 'rv', 'G', 'Grp', 'Gbp'])
    return t_f
    
    
# Estimate the Grp to H or K color for quasars given the redshift.
cDict = pickle.load(open('{}/quasar_color.pck'.format(modulepath), 'rb'))

def color_GrpH(z):
    '''
    Color of GRP-H.
    '''
    zm = cDict['z']
    color = cDict['G-H']
    return np.interp(z, zm, color)
    

def color_GrpK(z):
    '''
    Color of GRP-K.
    '''
    zm = cDict['z']
    color = cDict['G-K']
    return np.interp(z, zm, color)