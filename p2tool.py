import numpy as np
from datetime import datetime

now = datetime.now()
dt_string = now.strftime("%Y-%m-%dT%H:%M:%S")

all = ['check_item_exist', 'create_OB', 'ob_add_target', 'ob_add_description', 
       'template_add_acq_dual_wide', 'template_add_dual_obs_exp']


def check_item_exist(name, itemType, runContainerId, api):
    '''
    Check whether folder is in the container.
    
    Parameters
    ----------
    name : string
        The name of the item.
    itemType : string
        The type of the item, 'OB' or 'Folder'.
    runContainerId : int
        The ID of the container.
    api : p2api
    
    Returns
    -------
    it : dict 
        The information of the iterm if the item exists.  False, otherwise.
    '''
    items, _ = api.getItems(runContainerId)
    
    for it in items:
        if (it['itemType'] == itemType) & (it['name'] == name):
            return it
        
    return False


def create_OB(ob_name, folder_name, runContainerId, api):
    '''
    Create an OB or replace an OB in a folder.
    
    Parameters
    ----------
    ob_name : string
        The OB name.
    folder_name : string
        The folder name.
    runContainerId : int
        The ID of the container.
    api : p2api
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    folder_info = check_item_exist(folder_name, 'Folder', runContainerId, api)

    if folder_info:
        folderId = folder_info['containerId']
    else:
        folder, folderVersion = api.createFolder(runContainerId, folder_name)
        folderId = folder['containerId']
    
    ob_info = check_item_exist(ob_name, 'OB', folderId, api)
    if ob_info:
        obId = ob_info['obId']
        ob, obVersion = api.getOB(obId)
        api.deleteOB(obId, obVersion)
    ob, obVersion = api.createOB(folderId, ob_name)
    return ob, obVersion
    
    
def ob_add_target(ob, name, ra, dec, pma, pmd):
    '''
    Add the target information.
    
    Parameters
    ----------
    ob : dict
        The OB information.
    name : string
        The target name.
    ra : string
        The R.A. in hour angle, HH:MM:SS.SSS.
    dec : string
        The Decl. in degree, DD:MM:SS.SSS.
    pma : float
        The proper motion of R.A., in milliarcsec.
    pmd : float
        The proper motion of Decl., in milliarcsec.
    parallax : float
        The parallax, in milliarcsec.
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    ob['target']['name'] = name
    ob['target']['ra'] = ra
    ob['target']['dec'] = dec
    ob['target']['properMotionRa'] = np.round(pma * 1e-3, 6)
    ob['target']['properMotionDec'] = np.round(pmd * 1e-3, 6)
    return ob


def ob_add_description(ob, name, userComments):
    '''
    Parameters
    ----------
    ob : dict
        The OB information.
    name : string
        The observing description name.
    userComments : string
        The user comments.
    
    Returns
    -------
    ob : dict
        The OB information.
    '''
    ob['obsDescription']['name'] = name
    ob['obsDescription']['userComments'] = userComments
    return ob
    

def template_add_acq_dual_wide(api, obId, ft_name, ra_ft, dec_ft, pma_ft, pmd_ft, plx_ft, G_ft, H_ft, K_ft, 
                               sc_name, plx_sc, K_sc, ft_mode='AUTO', met_mode='ON', res='MED', pol='IN', 
                               ao_type='ADAPT_OPT', baseline=None, vltitype=None):
    '''
    Add GRAVITY_dual_wide_acq template.
    
    Parameters
    ----------
    obId : int
        The OB ID.
    ft_name : string
        The name of the FT target.
    ra_ft : float
        The right ascension of the ft source (degree or HH:MM:SS).
    dec_ft : float
        The declination of the ft source (degree or DD:MM:SS).
    pma_ft : float
        Proper motion in mas.
    pmd_ft : float
        Proper motion in mas.
    plx_ft : float
        Parallax in mas.
    G_ft : float
        The optical magnitude.
    H_ft : float
        The H band magnitude.
    K_ft : float
        The K band magnitude.
    sc_name : string
        The name of the SC target.
    plx_sc : float
        Parallax in mas.
    K_sc : float
        The K magnitude of the SC target.
    ft_mode : string (default: 'AUTO')
        The fringe tracker mode, 'AUTO', '1', '2', '7', '9'.
    met_mode : string (default: 'OFF')
        The metrology mode, 'ON', 'OFF', 'FAINT'.
    res : string (default: 'MED') 
        The spectral resolution.
    pol : string (default: 'IN')
        The polarization modes, 'IN' or 'OUT'.
    ao_type : string (default: 'ADAPT_OPT')
        The type of the AO, 'DEFAULT', 'AUTO_GUIDE', 'ADAPT_OPT', 'ADAPT_OPT_TCCD', 'IR_AO_OFFAXIS'.
        
    Returns
    -------
    acqTpl : dict
        Template information.
    acqTplVersion : string
        Template version.
    '''
    acqTpl, acqTplVersion = api.createTemplate(obId, 'GRAVITY_dual_wide_acq')
    
    plx_ft = np.round(plx_ft * 1e-3, 6)
    pma_ft = np.round(pma_ft * 1e-3, 6)
    pmd_ft = np.round(pmd_ft * 1e-3, 6)
    plx_sc = np.round(plx_sc * 1e-3, 6)
    
    if baseline is None:
        baseline = ['astrometric']
        
    if vltitype is None:
        vltitype = ['astrometry']
    
    acqParams = {
        'SEQ.FT.ROBJ.NAME': ft_name,
        'SEQ.FT.ROBJ.MAG': K_ft,
        'SEQ.INS.SOBJ.NAME': sc_name,
        'SEQ.INS.SOBJ.MAG': K_sc,
        'SEQ.FI.HMAG': H_ft,
        'SEQ.FT.ROBJ.ALPHA': ra_ft,
        'SEQ.FT.ROBJ.DELTA': dec_ft,
        'SEQ.FT.ROBJ.PARALLAX': plx_ft,
        'SEQ.FT.ROBJ.PMA': pma_ft,
        'SEQ.FT.ROBJ.PMD': pmd_ft,
        'SEQ.FT.MODE': ft_mode,
        'SEQ.MET.MODE': met_mode,
        'TEL.TARG.PARALLAX': plx_sc,
        'INS.SPEC.RES': res,
        'INS.FT.POL': pol,
        'INS.SPEC.POL': pol,
        'COU.AG.GSSOURCE': 'SETUPFILE',
        'COU.AG.ALPHA': ra_ft,
        'COU.AG.DELTA': dec_ft,
        'COU.GS.MAG': G_ft,
        'COU.AG.PMA': pma_ft,
        'COU.AG.PMD': pmd_ft,
        'COU.AG.TYPE': ao_type,
        'ISS.BASELINE': baseline,
        'ISS.VLTITYPE': vltitype
    }
    acqTpl, acqTplVersion = api.setTemplateParams(obId, acqTpl, acqParams, acqTplVersion)
    return acqTpl, acqTplVersion
    

def template_add_dual_obs_exp(api, obId, dit, ndit, hwpoff=[0], off_x=[0], off_y=[0], 
                              sky_x=2000, sky_y=2000, obsseq='O S'):
    '''
    Add the template of GRAVITY_dual_obs_exp.
    
    Parameters
    ----------
    obId : int
        The OB ID.
    dit : int
        DIT, 0.3, 1, 3, 10, 30, 100.
    ndit : int
        NDIT. Better to be multiple of 4 or 8.
    obsseq : string
        A number of "O" or "S".
        
    Returns
    -------
    scTpl : dict
        Template information.
    scTplVersion : string
        Template version.
    '''
    scTpl, scTplVersion = api.createTemplate(obId, 'GRAVITY_dual_obs_exp')
    scTpl, scTplVersion = api.setTemplateParams(obId, scTpl, {
        'DET2.DIT' : dit,
        'DET2.NDIT.OBJECT': ndit,
        'DET2.NDIT.SKY': ndit,
        'SEQ.HWPOFF': hwpoff,
        'SEQ.OBSSEQ': obsseq,
        'SEQ.RELOFF.X': off_x,
        'SEQ.RELOFF.Y': off_y,
        'SEQ.SKY.X': sky_x,
        'SEQ.SKY.Y': sky_y
    }, scTplVersion)
    return scTpl, scTplVersion