import numpy as np
import datetime
from .utils import load_obd_module

__all__ = ['sequencer', 'initiate_obd', 'add_acquisition_dual', 'add_exposure_dual', 'add_swap']


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


def add_acquisition_dual(obd_dict, tpl_name='', acq_dit='0.7', ft_pos={}, sc_pos={}, ft_name='',
                         ft_kmag='', ft_mode='AUTO', sc_name='', sc_kmag='', sobj_x='', sobj_y='',
                         acq_hmag='', res='MED', pol='IN', gssrc='SCIENCE', ao_pos={}, ao_mag='',
                         ao_type='DEFAULT', baseline='astrometric', vltitype='astrometry'):
    '''
    Add the acquisition template of dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_name : string (default: '')
        The template name.
    acq_dit : string (default: '0.7')
        The DIT of the acquisition camera.
    ft_pos : dict (default: {})
        The position of the FT source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
    sc_pos : dict (default: {})
        The position of the SC source.
            'ra' ('000000.000')
            'dec' ('000000.000')
            'pma' ('0.000')
            'pmd' ('0.000')
            'parallax' ('0.0')
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
    ft_parallax = ft_pos.get('parallax', '0.0')

    # SC coordinates
    sc_ra = sc_pos.get('ra', '000000.000')
    sc_dec = sc_pos.get('dec', '000000.000')
    sc_pma = sc_pos.get('pma', '0.000')
    sc_pmd = sc_pos.get('pmd', '0.000')
    sc_parallax = sc_pos.get('parallax', '0.0')

    # AO coordinates
    ao_ra = ao_pos.get('ra', '000000.000')
    ao_dec = ao_pos.get('dec', '000000.000')
    ao_pma = ao_pos.get('pma', '0.000')
    ao_pmd = ao_pos.get('pmd', '0.000')

    d = obd_dict['GRAVITY_dual_acq']
    d['TPL.NAME'] = tpl_name
    d['DET1.DIT'] = acq_dit
    d['TEL.TARG.ALPHA'] = ft_ra
    d['TEL.TARG.DELTA'] = ft_dec
    d['TEL.TARG.PMA'] = ft_pma
    d['TEL.TARG.PMD'] = ft_pmd
    d['TEL.TARG.PARALLAX'] = ft_parallax
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

    return obd_dict


def add_exposure_dual(obd_dict, tpl_num, tpl_name='', dit='1', ndit='8', hwpoff='0.0', obsseq='O S',
                      reloff_x='0.0', reloff_y='0.0', sky_x='2000', sky_y='2000'):
    '''
    Add the exposure of the dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_num : string
        The number of the template. It has to be unique for each template.
    tpl_name : string (default: '')
        The template name.
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
    
    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    key = '{0}_GRAVITY_dual_obs_exp'.format(tpl_num)
    if key in obd_dict:
        raise KeyError('There is already a "{}" in the obd_dict!'.format(key))
    
    obd_dict[key] = load_obd_module('GRAVITY_dual_obs_exp')
    obd_dict[key]['TPL.NAME'] = tpl_name
    obd_dict[key]['DET2.DIT'] = dit
    obd_dict[key]['DET2.NDIT.OBJECT'] = ndit
    obd_dict[key]['DET2.NDIT.SKY'] = ndit
    obd_dict[key]['SEQ.HWPOFF'] = hwpoff
    obd_dict[key]['SEQ.OBSSEQ'] = obsseq
    obd_dict[key]['SEQ.RELOFF.X'] = reloff_x
    obd_dict[key]['SEQ.RELOFF.Y'] = reloff_y
    obd_dict[key]['SEQ.SKY.X'] = sky_x
    obd_dict[key]['SEQ.SKY.Y'] = sky_y

    return obd_dict


def add_swap(obd_dict, tpl_num):
    '''
    Add the swap of the dual field observation.
    
    Parameters
    ----------
    obd_dict : string
        The OBD dict.
    tpl_num : string
        The number of the template. It has to be unique for each template.

    Returns
    -------
    obd_dict : dict
        The input OBD dict.
    '''
    key = '{0}_GRAVITY_dual_obs_swap'.format(tpl_num)
    if key in obd_dict:
        raise KeyError('There is already a "{}" in the obd_dict!'.format(key))
        
    obd_dict[key] = load_obd_module('GRAVITY_dual_obs_swap')
    return obd_dict


