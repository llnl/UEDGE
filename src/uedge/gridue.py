# Python routine for reading/writing HDF5-type grid files
# Routines are called from basis in favor of Fortran routines
from h5py import File

# TODO: Check that files exist, raise errors, flag for calling xerrab

def read_gridpars(fname=None):
    from uedge import bbb, com, grd
    dblxpts = ['dnull', 'snowflake15', 'snowflake45', 'snowflake75', 
        'dnXtarget', 'isoleg']
    geometry = com.geometry[0].decode('UTF-8').strip()
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()
    try:
        gridue = File(fname, 'r')
    except:
        gridue = File('{}.hdf5'.format(fname), 'r')
    _com = gridue['grid']['com']
    _grd = gridue['grid']['grd']
    com.nxm = _com['nxm'][()]
    com.nym = _com['nym'][()]
    if geometry in dblxpts:
        com.iysptrx1 = _com['iysptrx1'][()]
        com.iysptrx2 = _com['iysptrx2'][()]
        com.ixlb = _com['ixlb'][()]
        com.ixpt1 = _com['ixpt1'][()]
        com.ixmdp = _com['ixmdp'][()]
        com.ixpt2 = _com['ixpt2'][()]
        com.ixrb = _com['ixrb'][()]
        if geometry == 'dnXtarget':
            com.nxc = com.ixmdp[0]
    else:
        com.ixpt1[0] = _com['ixpt1'][()]
        com.ixpt2[0] = _com['ixpt2'][()]
        com.iysptrx1[0] = _com['iysptrx1'][()]
        com.iysptrx2[0] = com.iysptrx1[0]
        com.ixlb[0] = 0
        com.ixrb[0] = com.nxm
    try:
        com.simagxs = _com['simagxs'][()]
    except:
        pass
    try:
        com.sibdrys = _com['sibdrys'][()]
    except:
        pass
    gridue.close()
        
    

def read_gridue(fname=None):
    from uedge import bbb, com, grd
    from Forthon import gchange
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()

    if len(fname.split('/'))>3:
        prfname="/".join(["..."]+fname.split("/")[-2:])

    if com.iprint != 0:
        print(' Reading grid data from {}.hdf5'.format(prfname))
    read_gridpars(fname)
    gchange('RZ_grid_info')
    try:
        gridue = File(fname, 'r')
    except:
        gridue = File('{}.hdf5'.format(fname), 'r')
    _com = gridue['grid']['com']
    _grd = gridue['grid']['grd']
    com.rm = _com['rm'][()]
    com.zm = _com['zm'][()]
    com.psi = _com['psi'][()]
    com.br = _com['br'][()]
    com.bz = _com['bz'][()]
    com.bpol = _com['bpol'][()]
    com.bphi = _com['bphi'][()]
    com.b = _com['b'][()]
    com.runid = _com['runid'][()]
    try:
        com.nlim = _com['nlim'][()]
        gchange('Comflxgrd')
        com.xlim = _com['xlim'][()]
        com.ylim = _com['ylim'][()]
    except:
        pass
    try:
        grd.nplate1 = _grd['nplate1'][()]
        gchange('Mmod')
        grd.rplate1 = _grd['rplate1'][()]
        grd.zplate1 = _grd['zplate1'][()]
    except:
        pass
    try:
        grd.nplate2 = _grd['nplate2'][()]
        gchange('Mmod')
        grd.rplate2 = _grd['rplate2'][()]
        grd.zplate2 = _grd['zplate2'][()]
    except:
        pass
    gridue.close()
    if com.iprint != 0:
        print(' Grid data read successfully:')
        print('     file name:   {}.hdf5'.format(prfname))
        print('     run-ID:      {}'.format(com.runid[0].decode('UTF-8')))

def write_gridue(fname=None, runid=None):
    ''' Writes HDF5 grid file with name GridFileName '''
    from uedge import bbb, com, grd
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()
    if runid is None:
        runid = com.runid[0].decode('UTF-8').strip()
    if 'hdf5' not in fname.lower():
        fname = fname + '.hdf5'
    gridue = File(fname, 'w')
    gridue.require_group('grid')
    gue = gridue['grid']
    gue.require_group('com')
    gue.require_group('grd')
    _com = gue['com']
    _grd = gue['grd']

    _com.create_dataset('nxm', data=com.nxm)
    _com.create_dataset('nym', data=com.nym)
    _com.create_dataset('rm', data=com.rm)
    _com.create_dataset('zm', data=com.zm)
    _com.create_dataset('psi', data=com.psi)
    _com.create_dataset('br', data=com.br)
    _com.create_dataset('bz', data=com.bz)
    _com.create_dataset('bpol', data=com.bpol)
    _com.create_dataset('bphi', data=com.bphi)
    _com.create_dataset('b', data=com.b)
    _com.create_dataset('runid', data=runid)
    
    if com.geometry[0].decode('UTF-8').strip() == 'dnull':
        _com.create_dataset('ixpt1', data=com.ixpt1)
        _com.create_dataset('ixpt2', data=com.ixpt2)
        _com.create_dataset('iysptrx1', data=com.iysptrx1)
        _com.create_dataset('iysptrx2', data=com.iysptrx2)
        _com.create_dataset('ixlb', data=com.ixlb)
        _com.create_dataset('ixmdp', data=com.ixmdp)
        _com.create_dataset('ixrb', data=com.ixrb)
    else:
        _com.create_dataset('ixpt1', data=com.ixpt1[0])
        _com.create_dataset('ixpt2', data=com.ixpt2[0])
        _com.create_dataset('iysptrx1', data=com.iysptrx1[0])
        
    
    # Store extra data, such as limiter and plate data
    try:
        _com.create_dataset('simagxs', data=com.simagxs)
    except:
        pass
    try:
        _com.create_dataset('sibdrys', data=com.sibdrys)
    except:
        pass
    try:
        _com.create_dataset('nlim', data=com.nlim)
        _com.create_dataset('xlim', data=com.xlim)
        _com.create_dataset('ylim', data=com.ylim)
    except:
        pass
    try:
        _grd.create_dataset('nplate1', data=grd.nplate1)
        _grd.create_dataset('rplate1', data=grd.rplate1)
        _grd.create_dataset('zplate1', data=grd.zplate1)
    except:
        pass
    try:
        _grd.create_dataset('nplate2', data=grd.nplate1)
        _grd.create_dataset('rplate2', data=grd.rplate1)
        _grd.create_dataset('zplate2', data=grd.zplate1)
    except:
        pass
     

    gridue.close()
    if com.iprint != 0:
        print(' Wrote grid file successfully:')
        print('     file name:   {}'.format(fname.strip()))
        print('     run-ID:      {}'.format(runid.strip()))
