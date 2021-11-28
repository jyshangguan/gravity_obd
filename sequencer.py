import numpy as np
import datetime
from astropy.coordinates import SkyCoord
from astropy import units as u
from .utils import *

__all__ = ['gen_OBD_dual', 'gen_OBD_single', 
           'sequencer', 'initiate_obd', 'add_acquisition_dual', 'add_acquisition_single',
           'add_exposure_dual', 'add_exposure_single', 'add_swap']


def gen_OBD_dual(obd_name, ra_ft, dec_ft, pma_ft, pmd_ft, parallax_ft, radvel_ft, G_ft, K_ft, H_ft, 
                 ra_sc, dec_sc, pma_sc, pmd_sc, parallax_sc, radvel_sc, K_sc, 
                 acq_dit, dit, ndit, obsseq, obsid='00001', runID='Test 001', 
                 ft_name='s_ft', sc_name='s_sc', ft_mode='AUTO', res='MED', pol='IN', 
                 ao_type='ADAPT_OPT', baseline='astrometric', vltitype='astrometry'):
    '''
    Generate the OBD for the dual field observation.
    '''
    # Pos dict
    ft_pos = get_pos_dict(ra=ra_ft, dec=dec_ft, pma=pma_ft, pmd=pmd_ft, parallax=plx_ft, radvel=radvel_ft)
    sc_pos = get_pos_dict(ra=ra_sc, dec=dec_sc, pma=pma_sc, pmd=pmd_sc, parallax=plx_sc, radvel=radvel_sc)
    
    # offset
    c_ft = SkyCoord(ra_ft, dec_ft, frame='icrs', unit='deg')
    c_sc = SkyCoord(ra_sc, dec_sc, frame='icrs', unit='deg')
    sobj_x, sobj_y = sc_offset(c_sc, c_ft, pma_sc=pma_sc, pmd_sc=pmd_sc, plx_sc=plx_sc, radvel_sc=radvel_sc,
                               pma_ft=pma_ft, pmd_ft=pmd_ft, plx_ft=plx_ft, radvel_ft=radvel_ft)
    sobj_x = np.round(sobj_x.value, decimals=3)
    sobj_y = np.round(sobj_y.value, decimals=3)
    
    acq_kwargs = dict(acq_dit=acq_dit, ft_pos=ft_pos, sc_pos=sc_pos, ft_name=ft_name, ft_kmag=K_ft, ft_mode=ft_mode, 
                      sc_name=sc_name, sc_kmag=K_sc, sobj_x=sobj_x, sobj_y=sobj_y, acq_hmag=H_ft, res=res, pol=pol, 
                      ao_pos=None, ao_mag=G_ft, ao_type=ao_type, baseline=baseline, vltitype=vltitype)
   
    obsseq = obsseq.strip()
    obsseq = obsseq.replace(' ', '')
    tplParList = []
    for obs in obsseq:
        if obs not in ['S', 'O']:
            raise KeyError('Cannot recognize {} in the sequence!'.format(obs))
        tplParList.append(('exposure_dual', dict(acq_dit=acq_dit, dit=dit, ndit=ndit, obsseq=obs)))
    
    obd_dict = sequencer(obsid=obsid, runID=runID, acq_mode='dual', acq_kwargs=acq_kwargs, tplParList=tplParList)
    write_obd(obd_name, obd_dict)


def gen_OBD_single(obd_name, ra_sc, dec_sc, pma_sc, pmd_sc, parallax_sc, radvel_sc, K_sc, H_sc, G_sc, 
                   acq_dit, dit, ndit, obsseq, obsid='00001', runID='Test 001', sc_name='s_sc', 
                   ft_mode='AUTO', res='MED', pol='IN', ao_type='ADAPT_OPT', baseline='astrometric', 
                   vltitype='astrometry'):
    '''
    Generate the OBD for the single field observation.
    '''
    # Pos dict
    sc_pos = get_pos_dict(ra=ra_sc, dec=dec_sc, pma=pma_sc, pmd=pmd_sc, parallax=parallax_sc, radvel=radvel_sc)
    
    # offset
    c_sc = SkyCoord(ra_sc, dec_sc, frame='icrs', unit='deg')
    
    acq_kwargs = dict(acq_dit=acq_dit, sc_pos=sc_pos, ft_mode=ft_mode, sc_name=sc_name, 
                      sc_kmag=K_sc, acq_hmag=H_sc, res=res, pol=pol, ao_pos=None, 
                      ao_mag=G_sc, ao_type=ao_type, baseline=baseline, vltitype=vltitype)
   
    obsseq = obsseq.strip()
    obsseq = obsseq.replace(' ', '')
    tplParList = []
    for obs in obsseq:
        if obs not in ['S', 'O']:
            raise KeyError('Cannot recognize {} in the sequence!'.format(obs))
        tplParList.append(('exposure_dual', dict(acq_dit=acq_dit, dit=dit, ndit=ndit, obsseq=obs)))
    
    obd_dict = sequencer(obsid=obsid, runID=runID, acq_mode='single', acq_kwargs=acq_kwargs, tplParList=tplParList)
    write_obd(obd_name, obd_dict)
    

def sequencer(obsid=None, runID='wait...', now_str=None, acq_mode='dual', acq_kwargs={}, tplParList=[]):
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


def initiate_obd(obsid=None, runID='wait for the update!!!', now_str=None):
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


def add_acquisition_dual(obd_dict, acq_dit='0.7', ft_pos={}, sc_pos={}, ft_name='', ft_kmag='', 
                         ft_mode='AUTO', sc_name='', sc_kmag='', sobj_x='', sobj_y='', acq_hmag='', 
                         res='MED', pol='IN', ao_pos=None, ao_mag='', ao_type='DEFAULT', 
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
    acqTPLs = ['GRAVITY_dual_acq', 'GRAVITY_single_acq']
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
                 'TEL.TARG.RADVEL', 'SEQ.FT.ROBJ.NAME', 'SEQ.FT.ROBJ.MAG', 'SEQ.FT.MODE', 'SEQ.INS.SOBJ.NAME',
                 'SEQ.INS.SOBJ.MAG', 'SEQ.INS.SOBJ.X', 'SEQ.INS.SOBJ.Y', 'SEQ.FI.HMAG', 'INS.SPEC.RES', 'INS.FT.POL', 
                 'INS.SPEC.POL', 'COU.AG.GSSOURCE', 'COU.AG.ALPHA', 'COU.AG.DELTA', 'COU.GS.MAG', 'COU.AG.PMA', 
                 'COU.AG.PMD']
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

    
def add_acquisition_single(obd_dict, acq_dit='0.7', ft_mode='AUTO', sc_pos={}, sc_name='', sc_kmag='', 
                           acq_hmag='', res='MED', pol='IN', ao_pos=None, ao_mag='', 
                           ao_type='DEFAULT', baseline='astrometric', vltitype='astrometry', 
                           add_par_dict={}):
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
    add_par_dict : dict
        The additional keywords and values that are used.
        
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    acqTPLs = ['GRAVITY_dual_acq', 'GRAVITY_single_acq']
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

