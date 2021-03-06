import numpy as np
import datetime
from astropy.coordinates import SkyCoord
from astropy import units as u
from .utils import *

__all__ = ['gen_OBD_dual', 'gen_OBD_single', 'gen_OBD_wide',
           'sequencer', 'initiate_obd', 'add_acquisition_dual', 'add_acquisition_single',
           'add_exposure_dual', 'add_exposure_single', 'add_swap']

def gen_OBD_wide(obd_name, ra_ft, dec_ft, pma_ft, pmd_ft, parallax_ft, radvel_ft, G_ft, K_ft, H_ft, 
                 ra_sc, dec_sc, pma_sc, pmd_sc, parallax_sc, radvel_sc, K_sc, H_sc, dit, ndit, obsseq, 
                 acq_dit=0.7, sky_x=2000, sky_y=2000, ft_name='s_ft', sc_name='s_sc', ft_mode='AUTO', 
                 res='MED', pol='IN', met_mode='OFF', seq_align='T', ao_type='ADAPT_OPT', 
                 obsid='00001', runID='60.A-9102(I)', 
                 baseline='UTs', vltitype='snapshot'):
    '''
    Generate the OBD for the wide field observation.
    
    Parameters
    ----------
    ra_ft : float
        The right ascension of the ft source (degree or HH:MM:SS).
    dec_ft : float
        The declination of the ft source (degree or DD:MM:SS).
    pma_ft : float
        Proper motion in mas.
    pmd_ft : float
        Proper motion in mas.
    parallax_ft : float
        Parallax in mas.
    radvel_ft : float
        Radial velocity in km/s.
    G_ft : float
        The optical magnitude.
    K_ft : float
        The K band magnitude.
    H_ft : float
        The H band magnitude.
    ra_sc : float
        The right ascension of the ft source (degree or HH:MM:SS).
    dec_sc : float
        The declination of the ft source (degree or DD:MM:SS).
    pma_sc : float
        Proper motion in mas.
    pmd_sc : float
        Proper motion in mas.
    parallax_sc : float
        Parallax in mas.
    radvel_sc : float
        Radial velocity in km/s.
    K_sc : float
        The K magnitude of the SC target.
    H_sc : float
        The H magnitude of the SC target.
    dit : float or list
        The DIT of the science exposure.
    ndit : int or list
        The number of DITs of the science exposure.
    obsseq : list
        The list of exposures.
    acq_dit : float (default: 0.7)
        The DIT on the acquisition camera, units: second.
    sky_x : float (default: 2000)
        The offset for sky on RA, units: mas.
    sky_y : float (default: 2000)
        The offset for sky on DEC, units: mas.
    ft_name : string (default: 's_ft')
        The name of the FT target.
    sc_name : string (default: 's_sc') 
        The name of the SC target.
    ft_mode : string (default: 'AUTO')
        The fringe tracker mode, 'AUTO', '1', '2', '7', '9'.
    res : string (default: 'MED') 
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization modes, 'IN' or 'OUT'.
    met_mode : string (default: 'OFF')
        The metrology mode, 'ON', 'OFF', 'FAINT'.
    ao_type : string (default: 'ADAPT_OPT')
        The type of the AO, 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
    obsid : string
        The ID of the OB.
    runID : string
        The ID of the run.
    baseline : string (default: 'astrometric')
        The baseline types, 'small', 'medium', 'large', 'astrometric', 'UTs'.
    vltitype : string (default: 'astrometry')
        Not important, 'snapshot', 'imaging', 'time_series', 'astrometry'.
    
    Returns
    -------
    None
    '''
    assert isinstance(pma_ft, float)
    assert isinstance(pmd_ft, float)
    assert isinstance(parallax_ft, float)
    assert isinstance(radvel_ft, float)
    assert isinstance(pma_sc, float)
    assert isinstance(pmd_sc, float)
    assert isinstance(parallax_sc, float)
    assert isinstance(radvel_sc, float)
    assert isinstance(K_ft, float)
    assert isinstance(H_ft, float)
    assert isinstance(G_ft, float)
    assert isinstance(K_sc, float)
    assert isinstance(acq_dit, float)
    
    # Pos dict
    ft_pos = get_pos_J2000(ra=ra_ft, dec=dec_ft, pma=pma_ft, pmd=pmd_ft, parallax=parallax_ft, radvel=radvel_ft)
    sc_pos = get_pos_J2000(ra=ra_sc, dec=dec_sc, pma=pma_sc, pmd=pmd_sc, parallax=parallax_sc, radvel=radvel_sc)
    ao_pos = get_pos_J2000(ra=ra_ft, dec=dec_ft, pma=pma_ft, pmd=pmd_ft, parallax=parallax_ft, radvel=radvel_ft)
    
    # offset -- Checked with Taro, no need to calculate the offset
    c_ft = read_coordinate(ra_ft, dec_ft)
    c_sc = read_coordinate(ra_sc, dec_sc)
    sobj_x, sobj_y = sc_offset(c_sc, c_ft, pma_sc=pma_sc, pmd_sc=pmd_sc, plx_sc=parallax_sc, radvel_sc=radvel_sc,
                               pma_ft=pma_ft, pmd_ft=pmd_ft, plx_ft=parallax_ft, radvel_ft=radvel_ft)
    sobj_x = np.round(sobj_x.value, decimals=3)
    sobj_y = np.round(sobj_y.value, decimals=3)
    
    acq_kwargs = dict(acq_dit='{0:.1f}'.format(acq_dit), 
                      ft_pos=ft_pos, 
                      sc_pos=sc_pos, 
                      ft_name=ft_name, 
                      ft_kmag='{0:.1f}'.format(K_ft), 
                      ft_mode=ft_mode, 
                      sc_name=sc_name, 
                      sc_kmag='{0:.1f}'.format(K_sc), 
                      sc_hmag='{0:.1f}'.format(H_sc), 
                      sobj_x='{0:.3f}'.format(sobj_x), 
                      sobj_y='{0:.3f}'.format(sobj_y), 
                      acq_hmag='{0:.1f}'.format(H_ft), 
                      sky_x='{0:.1f}'.format(sky_x),
                      sky_y='{0:.1f}'.format(sky_y),
                      seq_align=seq_align,
                      res=res, 
                      pol=pol, 
                      ao_pos=ao_pos, 
                      ao_mag='{0:.1f}'.format(G_ft), 
                      ao_type=ao_type, 
                      met_mode=met_mode, 
                      baseline=baseline, 
                      vltitype=vltitype)
   
    assert isinstance(obsseq, list), 'Require obsseq to be a list!'
    if isinstance(ndit, int):
        ditList = [dit for ii in range(len(obsseq))]
        nditList = [ndit for ii in range(len(obsseq))]
    elif isinstance(ndit, list):
        ditList = dit
        nditList = ndit
    else:
        raise ValueError('The type of ndit ({}) is not correct!'.format(type(ndit)))
        
    tplParList = []
    for obs, dit, ndit in zip(obsseq, ditList, nditList):
        # Do not change the acq_dit for the exposure
        tplParList.append(('exposure_dual', 
                           dict(acq_dit='0.7', dit=dit, ndit=ndit, obsseq=obs, 
                                sky_x=sky_x, sky_y=sky_y)))
    
    obd_dict = sequencer(obsid=obsid, runID=runID, acq_mode='wide', acq_kwargs=acq_kwargs, tplParList=tplParList)
    write_obd(obd_name, obd_dict)
    
    
def gen_OBD_dual(obd_name, ra_ft, dec_ft, pma_ft, pmd_ft, parallax_ft, radvel_ft, G_ft, K_ft, H_ft, 
                 sobj_x, sobj_y, K_sc, dit, ndit, obsseq, acq_dit=0.7, reloff_x=0.0, reloff_y=0.0, 
                 sky_x=2000, sky_y=2000, ft_name='s_ft', sc_name='s_sc', ft_mode='AUTO', res='MED', 
                 pol='IN', met_mode='ON', ao_type='ADAPT_OPT', obsid='00001', runID='60.A-9102(I)', 
                 baseline='astrometric', vltitype='astrometry', gravity_mode='default', 
                 dit_on=None, ndit_on=None):
    '''
    Generate the OBD for the dual field observation.
    
    Parameters
    ----------
    ra_ft : float
        The right ascension of the ft source (degree or HH:MM:SS).
    dec_ft : float
        The declination of the ft source (degree or DD:MM:SS).
    pma_ft : float
        Proper motion in mas.
    pmd_ft : float
        Proper motion in mas.
    parallax_ft : float
        Parallax in mas.
    radvel_ft : float
        Radial velocity in km/s.
    G_ft : float
        The optical magnitude.
    K_ft : float
        The K band magnitude.
    H_ft : float
        The H band magnitude.
    sobj_x : float
        The RA offset of SC target from the FT target, units: mas.
    sobj_y : float
        The DEC offset of SC target from the FT target, units: mas.
    K_sc : float
        The K magnitude of the SC target.
    dit : float
        The DIT of the science exposure.
    ndit : int
        The number of DITs of the science exposure.
    obsseq : list
        The list of exposures.
    acq_dit : float (default: 0.7)
        The DIT on the acquisition camera, units: second.
    reloff_x : float (default: 0.0)
        The SC fiber offset, units: mas.
    reloff_y : float (default: 0.0)
        The SC fiber offset, units: mas.
    sky_x : float (default: 2000)
        The offset for sky on RA, units: mas.
    sky_y : float (default: 2000)
        The offset for sky on DEC, units: mas.
    ft_name : string (default: 's_ft')
        The name of the FT target.
    sc_name : string (default: 's_sc') 
        The name of the SC target.
    ft_mode : string (default: 'AUTO')
        The fringe tracker mode, 'AUTO', '1', '2', '7', '9'.
    res : string (default: 'MED') 
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization modes, 'IN' or 'OUT'.
    met_mode : string (default: 'OFF')
        The metrology mode, 'ON', 'OFF', 'FAINT'.
    ao_type : string (default: 'ADAPT_OPT')
        The type of the AO, 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
    obsid : string
        The ID of the OB.
    runID : string
        The ID of the run.
    baseline : string (default: 'astrometric')
        The baseline types, 'small', 'medium', 'large', 'astrometric', 'UTs'.
    vltitype : string (default: 'astrometry')
        Not important, 'snapshot', 'imaging', 'time_series', 'astrometry'.
    gravity_mode : string (default: 'default')
        The observation mode of GRAVITY, 'default' or 'exoGRAVITY'. The 'exoGRAVITY' mode is to 
        do dual-field on-axis observation like exoGRAVITY.
    dit_on (optional): float
        The DIT of science exposure on the FT star, only used when gravity_mode='exoGRAVITY'.
    ndit_on (optional): int
        The NDIT of science exposure on the FT star, only used when gravity_mode='exoGRAVITY'.
    
    Notes
    -----
    We need to provide the on-sky position of the targets for normal dual field mode.
    '''
    if not (isinstance(pma_ft, float) or isinstance(pma_ft, str)):
        pma_ft = 0
    if not (isinstance(pmd_ft, float) or isinstance(pmd_ft, str)):
        pmd_ft = 0
    if not (isinstance(parallax_ft, float) or isinstance(parallax_ft, str)):
        parallax_ft = 0
    if not (isinstance(radvel_ft, float) or isinstance(radvel_ft, str)):
        radvel_ft = 0
    if not (isinstance(K_ft, float) or isinstance(K_ft, str)):
        K_ft = 0
    if not (isinstance(H_ft, float) or isinstance(H_ft, str)):
        H_ft = 0
    if not (isinstance(G_ft, float) or isinstance(G_ft, str)):
        G_ft = 0
    if not (isinstance(K_sc, float) or isinstance(K_sc, str)):
        K_sc = 0
        
    # Pos dict
    #ft_pos = get_pos_current(ra=ra_ft, dec=dec_ft, pma=pma_ft, pmd=pmd_ft, parallax=parallax_ft, radvel=radvel_ft)
    ft_pos = get_pos_J2000(ra=ra_ft, dec=dec_ft, pma=pma_ft, pmd=pmd_ft, parallax=parallax_ft, radvel=radvel_ft)
    
    if gravity_mode == 'default':
        sobj_x_init = sobj_x
        sobj_y_init = sobj_y
    elif gravity_mode == 'exoGRAVITY':
        sep = np.sqrt(sobj_x**2 + sobj_y**2)
        if sep > 150:
            sobj_x_init = sobj_x * 1e-3
            sobj_y_init = sobj_y * 1e-3
        else:
            sobj_x_init = sobj_x * 1e-2
            sobj_y_init = sobj_y * 1e-2
    else:
        raise KeyError('Cannot recognize the gravity_mode ({})!'.format(gravity_mode))
    
    acq_kwargs = dict(acq_dit='{0:.1f}'.format(acq_dit), 
                      ft_pos=ft_pos, 
                      ft_name=ft_name, 
                      ft_kmag='{0:.1f}'.format(K_ft), 
                      ft_mode=ft_mode, 
                      sc_name=sc_name, 
                      sc_kmag='{0:.1f}'.format(K_sc), 
                      sobj_x='{0:.3f}'.format(sobj_x_init), 
                      sobj_y='{0:.3f}'.format(sobj_y_init), 
                      acq_hmag='{0:.1f}'.format(H_ft), 
                      sky_x='{0:.1f}'.format(sky_x),
                      sky_y='{0:.1f}'.format(sky_y),
                      res=res, 
                      pol=pol, 
                      met_mode=met_mode, 
                      ao_pos=None, 
                      ao_mag=G_ft, 
                      ao_type=ao_type, 
                      baseline=baseline, 
                      vltitype=vltitype)
    
    # Convert to string
    for k in ['sobj_x', 'sobj_y', 'ft_kmag', 'sc_kmag', 'acq_hmag', 'ao_mag']:
        v = acq_kwargs[k]
        if isinstance(v, float):
            acq_kwargs[k] = '{0:.3f}'.format(v)
    
    assert isinstance(obsseq, list), 'Require obsseq to be a list!'
    if isinstance(ndit, int):
        ditList = [dit for ii in range(len(obsseq))]
        nditList = [ndit for ii in range(len(obsseq))]
    elif isinstance(ndit, list):
        ditList = dit
        nditList = ndit
    else:
        raise ValueError('The type of ndit ({}) is not correct!'.format(type(ndit)))
    
    tplParList = []
    if gravity_mode == 'default':
        for obs, dit, ndit in zip(obsseq, ditList, nditList):
            if obs == 'swap':
                tplParList.append(('swap', dict()))
            else:
                tplParList.append(('exposure_dual', 
                                   dict(acq_dit='0.7', dit=dit, ndit=ndit, obsseq=obs, 
                                        reloff_x=reloff_x, reloff_y=reloff_y, sky_x=sky_x, 
                                        sky_y=sky_y)))
    else:
        assert dit_on is not None, 'Need to specify dit_on!'
        assert ndit_on is not None, 'Need to specify ndit_on!'
        for obs, dit, ndit in zip(obsseq, ditList, nditList):
            # On FT, O
            if obs == 'on:O':
                tplParList.append(('exposure_dual', 
                                   dict(acq_dit='0.7', dit=dit_on, ndit=ndit_on, obsseq='O', 
                                        reloff_x=-sobj_x_init, reloff_y=-sobj_y_init, sky_x=sky_x, 
                                        sky_y=sky_y)))
            # On FT, O S
            elif obs == 'on:OS':
                tplParList.append(('exposure_dual', 
                                   dict(acq_dit='0.7', dit=dit_on, ndit=ndit_on, obsseq='O S', 
                                        reloff_x=-sobj_x_init, reloff_y=-sobj_y_init, sky_x=sky_x, 
                                        sky_y=sky_y)))
            elif obs == 'swap':
                raise KeyError('Cannot do swap for dual-field on-axis!')
            else:
                tplParList.append(('exposure_dual', 
                                   dict(acq_dit='0.7', dit=dit, ndit=ndit, obsseq=obs, 
                                        reloff_x=(sobj_x-sobj_x_init), reloff_y=(sobj_y-sobj_y_init), 
                                        sky_x=sky_x, sky_y=sky_y)))
    
    obd_dict = sequencer(obsid=obsid, runID=runID, acq_mode='dual', acq_kwargs=acq_kwargs, tplParList=tplParList)
    write_obd(obd_name, obd_dict)


def gen_OBD_single(obd_name, ra_sc, dec_sc, pma_sc, pmd_sc, parallax_sc, radvel_sc, 
                   K_sc, H_sc, G_sc, acq_dit, dit, ndit, obsseq, sky_x=2000, sky_y=2000, 
                   obsid='00001', runID='60.A-9102(I)', sc_name='s_sc', ft_mode='AUTO', 
                   res='MED', pol='IN', met_mode='ON', ao_type='ADAPT_OPT', 
                   baseline='astrometric', vltitype='astrometry'):
    '''
    Generate the OBD for the single field observation.
    '''
    if not (isinstance(pma_sc, float) or isinstance(pma_sc, str)):
        pma_sc = 0
    if not (isinstance(pmd_sc, float) or isinstance(pmd_sc, str)):
        pmd_sc = 0
    if not (isinstance(parallax_sc, float) or isinstance(parallax_sc, str)):
        parallax_sc = 0
    if not (isinstance(radvel_sc, float) or isinstance(radvel_sc, str)):
        radvel_sc = 0
    if not (isinstance(K_sc, float) or isinstance(K_sc, str)):
        K_sc = 0
    if not (isinstance(H_sc, float) or isinstance(H_sc, str)):
        H_sc = 0
    if not (isinstance(G_sc, float) or isinstance(G_sc, str)):
        G_sc = 0

    # Pos dict
    #sc_pos = get_pos_current(ra=ra_sc, dec=dec_sc, pma=pma_sc, pmd=pmd_sc, parallax=parallax_sc, radvel=radvel_sc)
    sc_pos = get_pos_J2000(ra=ra_sc, dec=dec_sc, pma=pma_sc, pmd=pmd_sc, parallax=parallax_sc, radvel=radvel_sc)
    
    # offset
    c_sc = SkyCoord(ra_sc, dec_sc, frame='icrs', unit='deg')
    
    acq_kwargs = dict(acq_dit='{0:.1f}'.format(acq_dit), 
                      sc_pos=sc_pos, 
                      ft_mode=ft_mode, 
                      sc_name=sc_name, 
                      sc_kmag='{0:.1f}'.format(K_sc), 
                      acq_hmag='{0:.1f}'.format(H_sc), 
                      sky_x='{0:.1f}'.format(sky_x),
                      sky_y='{0:.1f}'.format(sky_y),
                      res=res, 
                      pol=pol, 
                      ao_pos=None, 
                      ao_mag='{0:.1f}'.format(G_sc), 
                      ao_type=ao_type, 
                      met_mode=met_mode, 
                      baseline=baseline, 
                      vltitype=vltitype)

    # Convert to string
    for k in ['sc_kmag', 'acq_hmag', 'ao_mag']:
        v = acq_kwargs[k]
        if isinstance(v, float):
            acq_kwargs[k] = '{0:.3f}'.format(v)
    
    assert isinstance(obsseq, list), 'Require obsseq to be a list!'
    if isinstance(ndit, int):
        ditList = [dit for ii in range(len(obsseq))]
        nditList = [ndit for ii in range(len(obsseq))]
    elif isinstance(ndit, list):
        ditList = dit
        nditList = ndit
    else:
        raise ValueError('The type of ndit ({}) is not correct!'.format(type(ndit)))
    
    tplParList = []
    for obs, dit, ndit in zip(obsseq, ditList, nditList):
        if obs not in ['S', 'O']:
            raise KeyError('Cannot recognize {} in the sequence!'.format(obs))
        tplParList.append(('exposure_single', 
                           dict(acq_dit='0.7', dit=dit, ndit=ndit, obsseq=obs, 
                                sky_x=sky_x, sky_y=sky_y)))

    obd_dict = sequencer(obsid=obsid, runID=runID, acq_mode='single', acq_kwargs=acq_kwargs, tplParList=tplParList)
    write_obd(obd_name, obd_dict)


def sequencer(obsid=None, runID='60.A-9102(I)', now_str=None, acq_mode='dual', acq_kwargs={}, tplParList=[]):
    '''
    The sequencer of the OBD.
    
    Parameters
    ----------
    obsid (optional) : string
        The OBS ID.
    runID : string (default: wait for the update!!!)
        The program ID.
    now_str (optional) : string 
        The current date and time.
    acq_mode : string (default: 'dual')
        The mode of the observation, 'single', 'dual', or 'wide'.
    acq_kwargs : dict (default: {})
        The keywords of the function to add the acquisition template.
    tplParList : list (default: [])
        The list of templates and their parameters, (template name, parameter dict).
        
    Returns
    -------
    obd_dict : dict
        The dict of the OBD information, ready to write.
    '''
    if now_str is None:
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%dT%H:%M:%S")
    
    obd_dict = initiate_obd(obsid=obsid, runID=runID, now_str=now_str)
    
    if acq_mode == 'dual':
        add_acquisition_dual(obd_dict, **acq_kwargs)
    elif acq_mode == 'single':
        add_acquisition_single(obd_dict, **acq_kwargs)
    elif acq_mode == 'wide':
        add_acquisition_wide(obd_dict, **acq_kwargs)
    else:
        raise KeyError('The acq_mode ({}) is not recognized!'.format(acq_mode))
    
    for (ii, (tpl, tpl_kwargs)) in enumerate(tplParList):
        if tpl == 'exposure_dual':
            add_exposure_dual(obd_dict, ii, **tpl_kwargs)
        elif tpl == 'exposure_single':
            add_exposure_single(obd_dict, ii, **tpl_kwargs)
        elif tpl == 'swap':
            add_swap(obd_dict, ii)
        else:
            raise KeyError('Cannot recognize the template ({0})!'.format(kk))
    return obd_dict


def initiate_obd(obsid=None, runID='60.A-9102(I)', now_str=None):
    '''
    Add the PAF and OBS modules.
    
    Parameters
    ----------
    obsid (optional) : string
        The OBS ID.
    runID : string (default: wait for the update!!!)
        The program ID.
    now_str (optional) : string 
        The current date and time.
    
    Returns
    -------
    obd_dict : dict
        The dict with PAF and OBS information.
    '''
    if now_str is None:
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%dT%H:%M:%S")

    obd_dict = {}
    obd_dict['PAF'] = load_obd_module('paf')
    obd_dict['PAF']['PAF.CRTE.DAYTIM'] = now_str
    obd_dict['PAF']['PAF.LCHG.DAYTIM'] = now_str

    obd_dict['OBS info'] = load_obd_module('obs')
    if obsid is None:
        obsid = '-1'
    obd_dict['OBS info']['OBS.ID'] = obsid
    obd_dict['OBS info']['OBS.PROG.ID'] = runID

    return obd_dict


def add_acquisition_wide(obd_dict, acq_dit='0.7', ft_pos={}, sc_pos={}, ft_name='', ft_kmag='', 
                         ft_mode='AUTO', sc_name='', sc_kmag='', sc_hmag='', sobj_x='', sobj_y='', 
                         acq_hmag='', res='MED', pol='IN', ao_pos=None, ao_mag='', ao_type='DEFAULT', 
                         sky_x='2000', sky_y='2000', met_mode='ON', seq_align='F', baseline='astrometric', 
                         vltitype='astrometry', add_par_dict={}):
    '''
    Add the acquisition template of dual field wide mode observation. The "TEL.TARG" 
    information is now for SC while "FT.TARG" information is for FT.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    ft_pos : dict (default: {})
        The position of the FT source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
            'radvel' ('0.0')
    sc_pos : dict (default: {})
        The position of the SC source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
            'radvel' ('0.0')
    ft_name : string (default: '')
        The name of the FT source.
    ft_kmag : string (default: '')
        The K magnitude of the FT source.
    ft_mode : string (default: 'AUTO')
        The fringe tracking mode (1, 2, 7, 9).
    sc_name : string (default: '')
        The name of the SC source.
    sc_kmag : string (default: '')
        The K magnitude of the SC source.
    sc_hmag : string (default: '')
        The H magnitude of the SC source.
    sobj_x : string (default: '')
        The RA offset of the SC source.
    sobj_y : string (default: '')
        The DEC offset of the SC source.
    acq_hmag : string (default: '')
        The H magnitude of the source on the acquisition camera.
    res : string (default: 'MED')
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization.  Both FT and SC.
    ao_pos : dict (default: {})
        The position of the AO source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
    ao_mag : string (default: '')
        The optical magnitude of the AO source.
    ao_type : string (default: 'DEFAULT')
        The type of the AO: 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    acqTPLs = ['GRAVITY_single_acq', 'GRAVITY_dual_acq', 'GRAVITY_dual_wide_acq']
    for tpl in acqTPLs:
        if tpl in obd_dict: 
            raise KeyError('There is already an acquisition template in obd_dict ({})!'.format(tpl))
        
    obd_dict['GRAVITY_dual_wide_acq'] = load_obd_module('GRAVITY_dual_wide_acq')

    # FT coordinates
    ft_ra = ft_pos.get('ra', '000000.000')
    ft_dec = ft_pos.get('dec', '000000.000')
    ft_pma = ft_pos.get('pma', '0.000')
    ft_pmd = ft_pos.get('pmd', '0.000')
    ft_plx = ft_pos.get('parallax', '0.0')
    ft_rv = ft_pos.get('radvel', '0.0')

    # SC coordinates
    sc_ra = sc_pos.get('ra', '000000.000')
    sc_dec = sc_pos.get('dec', '000000.000')
    sc_pma = sc_pos.get('pma', '0.000')
    sc_pmd = sc_pos.get('pmd', '0.000')
    sc_plx = sc_pos.get('parallax', '0.0')
    sc_rv = sc_pos.get('radvel', '0.0')

    # AO coordinates
    if ao_pos is None: 
        gssrc = 'SCIENCE'
        ao_ra = '000000.000'
        ao_dec = '000000.000'
        ao_pma = ft_pma
        ao_pmd = ft_pmd
    else:    
        gssrc = 'SETUPFILE'
        ao_ra = ao_pos.get('ra', '000000.000')
        ao_dec = ao_pos.get('dec', '000000.000')
        ao_pma = ao_pos.get('pma', '0.000')
        ao_pmd = ao_pos.get('pmd', '0.000')

    checkList = ['DET1.DIT', 'TEL.TARG.ALPHA', 'TEL.TARG.DELTA', 'TEL.TARG.PMA', 'TEL.TARG.PMD', 'TEL.TARG.PARALLAX', 
                 'TEL.TARG.RADVEL', 'SEQ.FT.ROBJ.ALPHA', 'SEQ.FT.ROBJ.DELTA', 'SEQ.FT.ROBJ.PMA', 'SEQ.FT.ROBJ.PMD', 
                 'SEQ.FT.ROBJ.PARALLAX', 'SEQ.FT.ROBJ.RADVEL', 'SEQ.FT.ROBJ.NAME', 'SEQ.FT.ROBJ.MAG', 'SEQ.FT.MODE', 
                 'SEQ.INS.SOBJ.NAME', 'SEQ.INS.SOBJ.MAG', 'SEQ.INS.SOBJ.HMAG', 'SEQ.INS.SOBJ.X', 'SEQ.INS.SOBJ.Y', 
                 'SEQ.FI.HMAG', 'SEQ.MET.MODE', 'INS.SPEC.RES', 'INS.FT.POL', 'INS.SPEC.POL', 'COU.AG.GSSOURCE', 
                 'COU.AG.ALPHA', 'COU.AG.DELTA', 'COU.GS.MAG', 'COU.AG.PMA', 'COU.AG.PMD']
    d = obd_dict['GRAVITY_dual_wide_acq']
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
    
    d['DET1.DIT'] = acq_dit
    d['TEL.TARG.ALPHA'] = sc_ra
    d['TEL.TARG.DELTA'] = sc_dec
    d['TEL.TARG.PMA'] = sc_pma
    d['TEL.TARG.PMD'] = sc_pmd
    d['TEL.TARG.PARALLAX'] = sc_plx
    d['TEL.TARG.RADVEL'] = sc_rv
    d['SEQ.FT.ROBJ.ALPHA'] = ft_ra
    d['SEQ.FT.ROBJ.DELTA'] = ft_dec
    d['SEQ.FT.ROBJ.PMA'] = ft_pma
    d['SEQ.FT.ROBJ.PMD'] = ft_pmd
    d['SEQ.FT.ROBJ.PARALLAX'] = ft_plx
    d['SEQ.FT.ROBJ.RADVEL'] = ft_rv
    
    d['SEQ.FT.ROBJ.NAME'] = ft_name
    d['SEQ.FT.ROBJ.MAG'] = ft_kmag
    d['SEQ.FT.MODE'] = ft_mode
    d['SEQ.INS.SOBJ.NAME'] = sc_name
    d['SEQ.INS.SOBJ.MAG'] = sc_kmag
    d['SEQ.INS.SOBJ.HMAG'] = sc_hmag
    d['SEQ.INS.SOBJ.X'] = sobj_x
    d['SEQ.INS.SOBJ.Y'] = sobj_y
    d['SEQ.FI.HMAG'] = acq_hmag
    d['SEQ.SKY.X'] = sky_x
    d['SEQ.SKY.Y'] = sky_y
    d['SEQ.ALIGN'] = seq_align
    d['SEQ.MET.MODE'] = met_mode
    d['INS.SPEC.RES'] = res
    d['INS.FT.POL'] = pol
    d['INS.SPEC.POL'] = pol
    
    d['COU.AG.GSSOURCE'] = gssrc
    d['COU.AG.ALPHA'] = ao_ra
    d['COU.AG.DELTA'] = ao_dec
    d['COU.GS.MAG'] = ao_mag
    d['COU.AG.PMA'] = ao_pma
    d['COU.AG.PMD'] = ao_pmd

    if ao_type not in ['DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS']:
        raise KeyError('Cannot recognize ao_type ({0})!'.format(ao_type))
    d['COU.AG.TYPE'] = ao_type

    if baseline not in ['small', 'medium', 'large', 'astrometric', 'UTs']:
        raise KeyError('Cannot recognize baseline ({0})!'.format(baseline))
    d['ISS.BASELINE'] = baseline

    if vltitype not in ['snapshot', 'imaging', 'time_series', 'astrometry']:
        raise KeyError('Cannot recognize vltitype ({0})!'.format(vltitype))
    d['ISS.VLTITYPE'] = vltitype
    
    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_dual_acq.obd!'.format(k))

    return obd_dict

    
def add_acquisition_dual(obd_dict, acq_dit='0.7', ft_pos={}, ft_name='', ft_kmag='', 
                         ft_mode='AUTO', sc_name='', sc_kmag='', sobj_x='', sobj_y='', 
                         acq_hmag='', sky_x='2000', sky_y='2000', res='MED', pol='IN', 
                         ao_pos=None, ao_mag='', ao_type='DEFAULT', met_mode='ON', 
                         baseline='astrometric', vltitype='astrometry', add_par_dict={}):
    '''
    Add the acquisition template of dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    ft_pos : dict (default: {})
        The position of the FT source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
            'radvel' ('0.0')
    ft_name : string (default: '')
        The name of the FT source.
    ft_kmag : string (default: '')
        The K magnitude of the FT source.
    ft_mode : string (default: 'AUTO')
        The fringe tracking mode (1, 2, 7, 9).
    sc_name : string (default: '')
        The name of the SC source.
    sc_kmag : string (default: '')
        The K magnitude of the SC source.
    sobj_x : string (default: '')
        The RA offset of the SC source.
    sobj_y : string (default: '')
        The DEC offset of the SC source.
    acq_hmag : string (default: '')
        The H magnitude of the source on the acquisition camera.
    res : string (default: 'MED')
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization.  Both FT and SC.
    ao_pos : dict (default: {})
        The position of the AO source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
    ao_mag : string (default: '')
        The optical magnitude of the AO source.
    ao_type : string (default: 'DEFAULT')
        The type of the AO: 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    acqTPLs = ['GRAVITY_single_acq', 'GRAVITY_dual_acq', 'GRAVITY_dual_wide_acq']
    for tpl in acqTPLs:
        if tpl in obd_dict: 
            raise KeyError('There is already an acquisition template in obd_dict ({})!'.format(tpl))
        
    obd_dict['GRAVITY_dual_acq'] = load_obd_module('GRAVITY_dual_acq')

    # FT coordinates
    ft_ra = ft_pos.get('ra', '000000.000')
    ft_dec = ft_pos.get('dec', '000000.000')
    ft_pma = ft_pos.get('pma', '0.000')
    ft_pmd = ft_pos.get('pmd', '0.000')
    ft_plx = ft_pos.get('parallax', '0.0')
    ft_rv = ft_pos.get('radvel', '0.0')

    # AO coordinates
    if ao_pos is None: 
        gssrc = 'SCIENCE'
        ao_ra = '000000.000'
        ao_dec = '000000.000'
        ao_pma = ft_pma
        ao_pmd = ft_pmd
    else:    
        gssrc = 'SETUPFILE'
        ao_ra = ao_pos.get('ra', '000000.000')
        ao_dec = ao_pos.get('dec', '000000.000')
        ao_pma = ao_pos.get('pma', '0.000')
        ao_pmd = ao_pos.get('pmd', '0.000')

    checkList = ['DET1.DIT', 'TEL.TARG.ALPHA', 'TEL.TARG.DELTA', 'TEL.TARG.PMA', 'TEL.TARG.PMD', 'TEL.TARG.PARALLAX', 
                 'TEL.TARG.RADVEL', 'SEQ.FT.ROBJ.NAME', 'SEQ.FT.ROBJ.MAG', 'SEQ.FT.MODE', 'SEQ.INS.SOBJ.NAME',
                 'SEQ.INS.SOBJ.MAG', 'SEQ.INS.SOBJ.X', 'SEQ.INS.SOBJ.Y', 'SEQ.FI.HMAG', 'SEQ.MET.MODE', 'INS.SPEC.RES', 
                 'INS.FT.POL', 'INS.SPEC.POL', 'COU.AG.GSSOURCE', 'COU.AG.ALPHA', 'COU.AG.DELTA', 'COU.GS.MAG', 
                 'COU.AG.PMA', 'COU.AG.PMD']
    d = obd_dict['GRAVITY_dual_acq']
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
    
    d['DET1.DIT'] = acq_dit
    d['TEL.TARG.ALPHA'] = ft_ra
    d['TEL.TARG.DELTA'] = ft_dec
    d['TEL.TARG.PMA'] = ft_pma
    d['TEL.TARG.PMD'] = ft_pmd
    d['TEL.TARG.PARALLAX'] = ft_plx
    d['TEL.TARG.RADVEL'] = ft_rv
    
    d['SEQ.FT.ROBJ.NAME'] = ft_name
    d['SEQ.FT.ROBJ.MAG'] = ft_kmag
    d['SEQ.FT.MODE'] = ft_mode
    d['SEQ.INS.SOBJ.NAME'] = sc_name
    d['SEQ.INS.SOBJ.MAG'] = sc_kmag
    d['SEQ.INS.SOBJ.X'] = sobj_x
    d['SEQ.INS.SOBJ.Y'] = sobj_y
    d['SEQ.FI.HMAG'] = acq_hmag
    d['SEQ.SKY.X'] = sky_x
    d['SEQ.SKY.Y'] = sky_y
    d['INS.SPEC.RES'] = res
    d['INS.FT.POL'] = pol
    d['INS.SPEC.POL'] = pol
    
    d['COU.AG.GSSOURCE'] = gssrc
    d['COU.AG.ALPHA'] = ao_ra
    d['COU.AG.DELTA'] = ao_dec
    d['COU.GS.MAG'] = ao_mag
    d['COU.AG.PMA'] = ao_pma
    d['COU.AG.PMD'] = ao_pmd

    if met_mode not in ['ON', 'FAINT']:
        raise KeyError('met_mode has to be "ON" or "FAINT"!')
    d['SEQ.MET.MODE'] = met_mode
    
    if ao_type not in ['DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS']:
        raise KeyError('Cannot recognize ao_type ({0})!'.format(ao_type))
    d['COU.AG.TYPE'] = ao_type

    if baseline not in ['small', 'medium', 'large', 'astrometric', 'UTs']:
        raise KeyError('Cannot recognize baseline ({0})!'.format(baseline))
    d['ISS.BASELINE'] = baseline

    if vltitype not in ['snapshot', 'imaging', 'time_series', 'astrometry']:
        raise KeyError('Cannot recognize vltitype ({0})!'.format(vltitype))
    d['ISS.VLTITYPE'] = vltitype
    
    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_dual_acq.obd!'.format(k))

    return obd_dict

    
def add_acquisition_single(obd_dict, acq_dit='0.7', ft_mode='AUTO', sc_pos={}, sc_name='', sc_kmag='', 
                           sky_x='2000', sky_y='2000', acq_hmag='', res='MED', pol='IN', ao_pos=None, 
                           ao_mag='', ao_type='DEFAULT', met_mode='ON', baseline='astrometric', 
                           vltitype='astrometry', add_par_dict={}):
    '''
    Add the acquisition template of dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    ft_mode : string (default: 'AUTO')
        The fringe tracking mode (1, 2, 7, 9).
    sc_pos : dict (default: {})
        The position of the SC source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
            'radvel' ('0.0')
    sc_name : string (default: '')
        The name of the SC source.
    sc_kmag : string (default: '')
        The K magnitude of the SC source.
    acq_hmag : string (default: '')
        The H magnitude of the source on the acquisition camera.
    res : string (default: 'MED')
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization.  Both FT and SC.
    ao_pos : dict (default: {})
        The position of the AO source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
    ao_mag : string (default: '')
        The optical magnitude of the AO source.
    ao_type : string (default: 'DEFAULT')
        The type of the AO: 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
    met_mode : string (default: 'ON')
        The metrology mode: 'ON', 'OFF', 'FAINT', 'WIDE' (not available).
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    acqTPLs = ['GRAVITY_single_acq', 'GRAVITY_dual_acq', 'GRAVITY_dual_wide_acq']
    for tpl in acqTPLs:
        if tpl in obd_dict: 
            raise KeyError('There is already an acquisition template in obd_dict ({})!'.format(tpl))
        
    obd_dict['GRAVITY_single_acq'] = load_obd_module('GRAVITY_single_acq')

    # SC coordinates
    sc_ra = sc_pos.get('ra', '000000.000')
    sc_dec = sc_pos.get('dec', '000000.000')
    sc_pma = sc_pos.get('pma', '0.000')
    sc_pmd = sc_pos.get('pmd', '0.000')
    sc_plx = sc_pos.get('parallax', '0.0')
    sc_rv = sc_pos.get('radvel', '0.0')

    # AO coordinates
    if ao_pos is None: 
        gssrc = 'SCIENCE'
        ao_ra = '000000.000'
        ao_dec = '000000.000'
        ao_pma = sc_pma
        ao_pmd = sc_pmd
    else:    
        gssrc = 'SETUPFILE'
        ao_ra = ao_pos.get('ra', '000000.000')
        ao_dec = ao_pos.get('dec', '000000.000')
        ao_pma = ao_pos.get('pma', '0.000')
        ao_pmd = ao_pos.get('pmd', '0.000')

    checkList = ['DET1.DIT', 'TEL.TARG.ALPHA', 'TEL.TARG.DELTA', 'TEL.TARG.PMA', 'TEL.TARG.PMD', 
                 'TEL.TARG.PARALLAX', 'TEL.TARG.RADVEL', 'SEQ.FT.MODE', 'SEQ.INS.SOBJ.NAME', 
                 'SEQ.INS.SOBJ.MAG', 'SEQ.FI.HMAG', 'INS.SPEC.RES', 'INS.FT.POL', 'INS.SPEC.POL', 
                 'COU.AG.GSSOURCE',  'COU.AG.ALPHA', 'COU.AG.DELTA', 'COU.GS.MAG', 'COU.AG.PMA', 
                 'COU.AG.PMD']
    d = obd_dict['GRAVITY_single_acq']
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
    
    d['DET1.DIT'] = acq_dit
    d['TEL.TARG.ALPHA'] = sc_ra
    d['TEL.TARG.DELTA'] = sc_dec
    d['TEL.TARG.PMA'] = sc_pma
    d['TEL.TARG.PMD'] = sc_pmd
    d['TEL.TARG.PARALLAX'] = sc_plx
    d['TEL.TARG.RADVEL'] = sc_rv
    
    d['SEQ.FT.MODE'] = ft_mode
    d['SEQ.INS.SOBJ.NAME'] = sc_name
    d['SEQ.INS.SOBJ.MAG'] = sc_kmag
    d['SEQ.FI.HMAG'] = acq_hmag
    d['SEQ.SKY.X'] = sky_x
    d['SEQ.SKY.Y'] = sky_y
    d['SEQ.MET.MODE'] = met_mode
    d['INS.SPEC.RES'] = res
    d['INS.FT.POL'] = pol
    d['INS.SPEC.POL'] = pol
    
    d['COU.AG.GSSOURCE'] = gssrc
    d['COU.AG.ALPHA'] = ao_ra
    d['COU.AG.DELTA'] = ao_dec
    d['COU.GS.MAG'] = ao_mag
    d['COU.AG.PMA'] = ao_pma
    d['COU.AG.PMD'] = ao_pmd

    if ao_type not in ['DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS']:
        raise KeyError('Cannot recognize ao_type ({0})!'.format(ao_type))
    d['COU.AG.TYPE'] = ao_type

    if baseline not in ['small', 'medium', 'large', 'astrometric', 'UTs']:
        raise KeyError('Cannot recognize baseline ({0})!'.format(baseline))
    d['ISS.BASELINE'] = baseline

    if vltitype not in ['snapshot', 'imaging', 'time_series', 'astrometry']:
        raise KeyError('Cannot recognize vltitype ({0})!'.format(vltitype))
    d['ISS.VLTITYPE'] = vltitype
    
    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_dual_acq.obd!'.format(k))

    return obd_dict

    
def add_exposure_dual(obd_dict, tpl_num, acq_dit='0.7', dit='0.3', ndit='32', hwpoff='0.0', obsseq='O S',
                      reloff_x='0.0', reloff_y='0.0', sky_x='2000', sky_y='2000', 
                      add_par_dict={}):
    '''
    Add the exposure of the dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_num : string
        The number of the template. It has to be unique for each template.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    dit : string (default: '1')
        The DIT of the science exposure.
    ndit : string (default: '8')
        The NDIT of the science exposure.
    obsseq : string (default: 'O S')
        The observation sequence.
    reloff_x : string (default: '0.0')
        The RA offset of SC fiber.
    reloff_y : string (default: '0.0')
        The DEC offset of SC fiber.
    sky_x : string (default: '2000')
        The RA offset of the sky.
    sky_y : string (default: '2000')
        The DEC offset of the sky.
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    key = '{0}_GRAVITY_dual_obs_exp'.format(tpl_num)
    if key in obd_dict:
        raise KeyError('There is already a "{}" in the obd_dict!'.format(key))
    
    obd_dict[key] = load_obd_module('GRAVITY_dual_obs_exp')
    
    checkList = ['DET1.DIT', 'DET2.DIT', 'DET2.NDIT.OBJECT', 'DET2.NDIT.SKY', 'SEQ.HWPOFF',
                 'SEQ.OBSSEQ', 'SEQ.RELOFF.X', 'SEQ.RELOFF.Y', 'SEQ.SKY.X', 'SEQ.SKY.Y']
    d = obd_dict[key]
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
    
    d['DET1.DIT'] = acq_dit
    d['DET2.DIT'] = dit
    d['DET2.NDIT.OBJECT'] = ndit
    d['DET2.NDIT.SKY'] = ndit
    d['SEQ.HWPOFF'] = hwpoff
    d['SEQ.OBSSEQ'] = obsseq
    d['SEQ.RELOFF.X'] = reloff_x
    d['SEQ.RELOFF.Y'] = reloff_y
    d['SEQ.SKY.X'] = sky_x
    d['SEQ.SKY.Y'] = sky_y

    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_dual_obs_exp.obd!'.format(k))

    return obd_dict


def add_exposure_single(obd_dict, tpl_num, acq_dit='0.7', dit='0.3', ndit='32', hwpoff='0.0', obsseq='O S',
                        sky_x='2000', sky_y='2000', add_par_dict={}):
    '''
    Add the exposure of the single field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_num : string
        The number of the template. It has to be unique for each template.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    dit : string (default: '1')
        The DIT of the science exposure.
    ndit : string (default: '8')
        The NDIT of the science exposure.
    obsseq : string (default: 'O S')
        The observation sequence.
    sky_x : string (default: '2000')
        The RA offset of the sky.
    sky_y : string (default: '2000')
        The DEC offset of the sky.
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    key = '{0}_GRAVITY_single_obs_exp'.format(tpl_num)
    if key in obd_dict:
        raise KeyError('There is already a "{}" in the obd_dict!'.format(key))
    
    obd_dict[key] = load_obd_module('GRAVITY_single_obs_exp')
    
    checkList = ['DET1.DIT', 'DET2.DIT', 'DET2.NDIT.OBJECT', 'DET2.NDIT.SKY', 'SEQ.HWPOFF',
                 'SEQ.OBSSEQ', 'SEQ.SKY.X', 'SEQ.SKY.Y']
    d = obd_dict[key]
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
    
    d['DET1.DIT'] = acq_dit
    d['DET2.DIT'] = dit
    d['DET2.NDIT.OBJECT'] = ndit
    d['DET2.NDIT.SKY'] = ndit
    d['SEQ.HWPOFF'] = hwpoff
    d['SEQ.OBSSEQ'] = obsseq
    d['SEQ.SKY.X'] = sky_x
    d['SEQ.SKY.Y'] = sky_y

    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_single_obs_exp.obd!'.format(k))

    return obd_dict


def add_swap(obd_dict, tpl_num, ft_mode='AUTO', add_par_dict={}):
    '''
    Add the swap of the dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_num : string
        The number of the template. It has to be unique for each template.
    ft_mode : string (default: 'AUTO')
        The fringe tracking mode (1, 2, 7, 9).
    add_par_dict : dict
        The additional keywords and values that are used.

    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    key = '{0}_GRAVITY_dual_obs_swap'.format(tpl_num)
    if key in obd_dict:
        raise KeyError('There is already a "{}" in the obd_dict!'.format(key))
        
    obd_dict[key] = load_obd_module('GRAVITY_dual_obs_swap')
    
    checkList = ['SEQ.FT.MODE']
    d = obd_dict[key]
    for k in checkList:
        if k not in d:
            raise KeyError('Cannot find {}!'.format(k))
            
    d['SEQ.FT.MODE'] = ft_mode
    
    for k in add_par_dict:
        if k in checkList:
            raise KeyError('The parameter ({}) is duplicated!'.format(k))
            
        if k in d:
            d[k] = add_par_dict[k]
        else:
            raise KeyError('The key ({}) is not found in GRAVITY_dual_obs_swap.obd!'.format(k))
            
    return obd_dict

