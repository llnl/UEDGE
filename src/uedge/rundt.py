# Holm10 Nov 5 2019, based on rdcontdt.py
# 191121 -  Created hdf5-routines to read and save time-dependent data
#           Writes and reads dictionary with multi-dimensional arrays
#           containing all restore-parameters.
# 230210 - Updated old rundt function to RunData class, modularizing
#          all functionalities. Added a comprehensive diagnostics suite
#          plotting fnrm evolution as function of exmains(), plasma 
#          time, and wall-clock time. Still need to test time-slicing
#   `      procedures
# 230327 - Removed old routine, created wrapper function rundt for 
#          object. Renamed Object to UeRun.
# 230522 - Fixed bug associated with itroub, improved itroub visualization

from matplotlib.pyplot import ion
ion()

class UeRun():
    ''' Class containing information on run '''
    def __init__(self, *args, n_stor = False, **kwargs):
        from numpy import array        
        from uedge import bbb, com
        # TODO: Add restore/recover from timeslice
        # TODO: Add plot timeslice directly
        # NOTE: No -> Utilize direct I/O from file instead
        self.numvar = bbb.numvar
        try:
            self.nx
        except:
            self.nx = com.nx
        try:
            self.ny
        except:
            self.ny = com.ny
        self.ixpt1 = com.ixpt1[0]
        self.ixpt2 = com.ixpt2[0]
        self.iysptrx = com.iysptrx
        self.equationkey = array([b'te', b'ti', b'phi', b'up', b'ni', b'ng', 
            b'tg'])

        self.setupvars = {}
        for var in [ 'numvar', 'isteon', 'istion', 'isupon', 'isphion', 
            'isupgon', 'isngon', 'istgon', 'ishymol', 'nisp', 'ngsp', 'nhsp', 
            'nhgsp', 'nzsp', 'b0', 'ncore', 'pcoree', 'pcorei', 'nx', 'ny',
            'iysptrx', 'ixpt1', 'ixpt2', 'neq'
            ]:
            try:
                self.setupvars[var] = getattr(bbb, var)
            except:
                self.setupvars[var] = getattr(com, var)
        self.setupvars['ixpt1'] = self.setupvars['ixpt1'][0]
        self.setupvars['ixpt2'] = self.setupvars['ixpt2'][0]        
        super().__init__(*args, **kwargs)

    def itroub(self):
        ''' Function that displays information on the problematic equation '''
        from numpy import mod, argmax, where, array, argmin
        from uedge import bbb
        from copy import deepcopy

        self.equations = [bbb.idxte, bbb.idxti, bbb.idxphi, 
            bbb.idxu, bbb.idxn, bbb.idxg, bbb.idxtg]
        equationsdescription = [ 'Electron energy', 'Ion energy', 'Potential',
            'Ion momentum', 'Ion density', 'Gas density', 'Gas temperature']
        # Assert required dict is present
        try:
            self.classvars
        except:
            self.classvars = {}
            for cvar in ['itrouble', 'troubleeq', 'internaleq', 
                    'internalspecies', 'troubleindex', 'dtrealfail',
                    'ylfail', 'yldotsfscalfail']:
                self.classvars[cvar] = []
        # Find the fortran index of the troublemaking equation
        self.classvars['itrouble'].append(deepcopy(argmax(abs(bbb.yldot*\
            bbb.sfscal)[:bbb.neq])+1))
        print("** Fortran index of trouble making equation is:\n{}".format(\
            self.classvars['itrouble'][-1]))
        # Print equation information
        print("** Number of equations solved per cell:\n    numvar = {}\n"\
            .format(bbb.numvar))

        self.classvars['troubleeq'].append(mod(self.classvars['itrouble'][-1]-1, bbb.numvar)+1)

        species = ''
        self.classvars['internaleq'].append([abs(x - self.classvars['itrouble'][-1]).min() for x in \
            self.equations].index(0))
        if self.equations[self.classvars['internaleq'][-1]].ndim == 3:
            self.classvars['internalspecies'].append( where(\
                self.equations[self.classvars['internaleq'][-1]] == self.classvars['itrouble'][-1])\
                [-1][0] + 1)
            species = ' of species {}'.format(self.classvars['internalspecies'][-1])
        else:
            self.classvars['internalspecies'].append(0)


        print('** Troublemaker equation is:\n{} equation{}: iv_t={}\n'\
            .format(equationsdescription[self.classvars['internaleq'][-1]], species, 
            self.classvars['troubleeq'][-1]))

        # Display additional information about troublemaker cell
        self.classvars['troubleindex'].append(deepcopy(bbb.igyl[self.classvars['itrouble'][-1]-1,]))
        self.classvars['dtrealfail'].append(deepcopy(bbb.dtreal))
        self.classvars['ylfail'].append(deepcopy(bbb.yl[self.classvars['itrouble'][-1]-1]))
        self.classvars['yldotsfscalfail'].append(deepcopy((bbb.yldot*bbb.sfscal)\
            [self.classvars['itrouble'][-1]-1]))
        print('** Troublemaker cell (ix,iy) is:\n' + \
            '{}\n'.format(self.classvars['troubleindex'][-1]))
        print('** Timestep for troublemaker equation:\n' + \
            '{:.4e}\n'.format(self.classvars['dtrealfail'][-1]))
        print('** yl for troublemaker equation:\n' + \
            '{:.4e}\n'.format(self.classvars['ylfail'][-1]))
        print('** yldot*sfscal for troublemaker equation:\n' + \
            '{:.4e}\n'.format(self.classvars['yldotsfscalfail'][-1]))
    
    def savesuccess(self, savename, fnrm=None):
        from time import time
        from uedge import bbb
        from copy import deepcopy

        try:
            self.classvars['time'].append(time())
        except:
            pass
        if 'fnorm' in self.classvars.keys():
            if fnrm is None:
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                self.classvars['fnorm'].append(deepcopy((sum((bbb.yldot[:bbb.neq]*\
                    bbb.sfscal[:bbb.neq])**2))**0.5))
            else:
                self.classvars['fnorm'].append(fnrm)
        for var in ['nfe', 'dt_tot', 'dtreal']:
            if var in self.classvars:
                self.classvars[var].append(deepcopy(getattr(bbb, var)))
        self.save_intermediate(savename)

    def store_timeslice(self):
        from copy import deepcopy
        from uedge import bbb

        for var in ['ni', 'ng', 'up', 'te', 'ti', 'tg', 'phi', 'dt_tot']:
            self.classvars['slice_{}'.format(var)].append(deepcopy(getattr(bbb, var)))

            
    def save_intermediate(self, savename):
        from uedge.hdf5 import hdf5_save
        from os.path import exists
        from uedge import bbb, com
        from h5py import File

        for var in [ 'isteon', 'istion', 'isupon', 'isphion', 'isupgon',
            'isngon', 'istgon', 'ishymol']:
            self.__setattr__(var, bbb.__getattribute__(var))
        for var in [ 'nisp', 'ngsp', 'nhsp', 'nhgsp', 'nzsp']:
            self.__setattr__(var, com.__getattribute__(var))
        if not exists('/'.join(savename.split('/')[:-1])):
            print('Folder {} not found, saving output to cwd...'\
                .format(savename.split('/')[0]))
            hdf5_save(savename)
        try:
            self.save(savename)#.replace('.hdf5', '_UeCase.hdf5'))
        except:
            hdf5_save(savename)

        with File(savename, 'r+') as file:
            file.require_group('convergence')
            group = file['convergence']
            group.create_dataset('t_start', data=self.tstart)
            group.create_dataset('equationkey', data=self.equationkey)
            for var, value in self.setupvars.items():
                group.create_dataset(var, data=value)
            for var, value in self.classsetup.items():
                group.create_dataset(var, data=value)
            for var in self.classvars:
                group.create_dataset(var, data=self.classvars[var])
        print('Intermediate solution written to {}'.format(savename))

    def convergenceanalysis(self, savefname, fig=None,
        xaxis = 'exmain', logx = False, color='k', label=None,
        ylim = (None, None)):
        from h5py import File
        from matplotlib.pyplot import subplots
        from os.path import exists
        from datetime import timedelta
        from matplotlib.ticker import FuncFormatter
        from numpy import cumsum, ones
        if fig is None:
            f, ax = subplots(1, 3, figsize=(15, 5))
        else:
            ax = fig.get_axes()
            if len(ax) < 3:
                print('Three subplots required for plots! Aborting...')
                return
            f = fig

        if not exists(savefname):
            print('File {} not found. Aborting!'.format(savefname))
            return

        with File(savefname, 'r') as file:
            data = file['convergence']
            try:
                data = file['convergence']
            except:
                print('Convergence data not found in {}. Aborting!'.format(\
                    savefname))
                return
            
            if xaxis == 'exmain':
                xlabel = 'Exmain calls'
                xones = ones(data['ii2'][()].shape)
                x = cumsum(xones)
            elif xaxis == 'nfe':
                xlabel = 'nfe internal iterations'
                x = cumsum(data['nfe'][()][:, 0, 0])
            elif xaxis == 'time':
                xlabel = 'Total wall-clock time [HH:MM]'
                x = [timedelta(t - data['t_start'][()]) for t in data['time'][()]]
                x = data['time'][()] - data['t_start'][()]
                
            if logx is True:
                ax[0].loglog(x, data['fnorm'][()], '-', color=color, label=label)
                ax[1].loglog(data['dt_tot'][()], data['fnorm'][()], '-', 
                    color=color, label=label)
                ax[2].loglog(x, data['dtreal'][()], '-', color=color, label=label)
            else:
                ax[0].semilogy(x, data['fnorm'][()], '-', color=color, label=label)
                ax[1].semilogy(data['dt_tot'][()], data['fnorm'][()], '-', 
                    color=color, label=label)
                ax[2].semilogy(x, data['dtreal'][()], '-', color=color, 
                    label=label)
            ax[1].set_title('Total exmain evaluations: {}'.format\
                (len(data['dtreal'][()])))

        ax[0].set_xlabel(xlabel)
        ax[1].set_xlabel('Accumulated plasma simualtion time [s]')
        ax[2].set_xlabel(xlabel)
        ax[0].set_ylabel('Initial fnorm')
        ax[1].set_ylabel('Initial fnorm')
        ax[2].set_ylabel('Time-step (dtreal) [s]')
        ax[0].set_ylim(ylim)
        ax[1].set_ylim(ylim)
        if xaxis == 'time':
            ax[0].xaxis.set_major_formatter(FuncFormatter(lambda t, pos :\
                    str(timedelta(seconds=t))[:-3]))
            ax[2].xaxis.set_major_formatter(FuncFormatter(lambda t, pos :\
                    str(timedelta(seconds=t))[:-3]))
        if label is not None:
            ax[0].legend()

        return f


    def failureanalysis(self, savefname, equation=None, N=slice(None), geometry=False):
        from h5py import File
        from os.path import exists
        from matplotlib.pyplot import subplots
        from numpy import histogram, zeros
        from matplotlib.collections import PolyCollection
        
        # TODO: add option for plotting offending cell on geometric grid
        if geometry is True:
            f, ax = subplots(1,2, figsize=(12, 8), 
                gridspec_kw={'width_ratios': [1, 1.2]})
            f.subplots_adjust(top=0.98)
        else:
            f, ax = subplots(2,1, figsize=(10, 7))

        if not exists(savefname):
            print('File {} not found. Aborting!'.format(savefname))
            return
        
        with File(savefname, 'r') as file:
            try:
                data = file['convergence']
            except:
                print('Convergence data not found in {}. Aborting!'.format(\
                    savefname))
                return

            if equation is not None:
                iequation = [x.decode('UTF-8') for x in data['equationkey']]\
                    .index(equation)
            # Bin the equation errors
            nspecies = 1/(data['nisp'][()] + 1)
            nbins = 7*data['nisp'][()]
            counts, bins = histogram((data['internaleq'][N]+\
                data['internalspecies'][N]*nspecies)-0.5, bins=nbins, 
                range=(-0.5,6.5))
            h, e = histogram(data['internaleq'][N] - 0.5, bins=7, 
                range=(-0.5,6.5))
            ax[0].bar([x for x in range(7)], h, width=1, edgecolor='k',
                color=(0, 87/255, 183/255))
            ax[0].hist(bins[3*data['nisp'][()]:-1], bins[3*data['nisp'][()]:], 
                weights=counts[3*data['nisp'][()]:], edgecolor='k', 
                color=(255/255, 215/255, 0))
            ax[0].set_xticks(range(7))
            ax[0].set_xticklabels([x.decode('UTF-8') for x in \
                data['equationkey'][()]])
            ax[0].grid(linestyle=':', linewidth=0.5, axis='y')
            ax[0].set_xlim((-0.5,6.5))
            ax[0].set_ylabel('Counts')
            for i in range(7):
                ax[0].axvline(i-0.5, linewidth=1, color='k')

            # Visualize error locations
            nx = data['nx'][()]        
            ny = data['ny'][()]        
            ixpt1 = data['ixpt1'][()]
            ixpt2 = data['ixpt2'][()]
            iysptrx = data['iysptrx'][()]
            frequency = zeros((nx+2, ny+2))

            cells = []
            if geometry is True:
                try:
                    rm = file['grid/com/rm'][()] 
                    zm = file['grid/com/zm'][()] 
                except:
                    raise KeyError('Grid data not found in {}'.format(\
                        savefname))
                for i in range(nx+2):
                    for j in range(ny+2):
                        cell = []
                        for k in [1, 2, 4, 3]:
                            cell.append((rm[i, j, k], zm[i, j, k]))
                        cells.append(cell)
            else:
                for i in range(nx+2):
                    for j in range(ny+2):
                        cells.append([[i-.5, j-.5], [i+.5, j-.5], 
                            [i+.5, j+.5], [i-.5, j+.5]])

            polys = PolyCollection(cells, edgecolors='k', 
                linewidth=0.5-0.45*geometry, linestyle=':')
                
            for i in range(len(data['itrouble'][N])):
                coord = data['troubleindex'][N][i]
                if equation is None:
                    frequency[coord[0], coord[1]] += 1
                elif iequation == data['internaleq'][N][i]:
                    frequency[coord[0], coord[1]] += 1

        if geometry is False:
            polys.set_cmap('binary')
        else:
            polys.set_cmap('binary')
        polys.set_array(frequency.reshape(((nx+2)*(ny+2),)))

        cbar = f.colorbar(polys, ax=ax[1])
        cbar.ax.set_ylabel('N trouble'+' for {}'.format(equation)*\
            (equation is not None), va='bottom', labelpad=20)


        ax[1].set_xlabel('Poloidal index')
        ax[1].set_ylabel('Radial index')
        ax[1].add_collection(polys)
        if geometry is False:
            ax[1].plot([.5, nx+.5, nx+.5, .5, .5], [.5, .5, ny+.5, ny+.5, .5], 
                'k-', linewidth=1)
            ax[1].plot([.5, nx+.5],[iysptrx+.5, iysptrx+.5], 'k-', 
                linewidth=1)
            ax[1].plot([ixpt1+.5, ixpt1+.5], [.5, iysptrx+.5], 'k-', 
                linewidth=1)
            ax[1].plot([ixpt2+.5, ixpt2+.5], [.5, iysptrx+.5], 'k-', 
                linewidth=1)
        else:
            ax[1].set_aspect('equal')
            ax[1].set_ylim((0.95*zm.min(), 1.05*zm.max()))
            ax[1].set_xlim((0.95*rm.min(), 1.05*rm.max()))

        return f


    def restorevalues(self):
        ''' Restores the original UEDGE values '''
        from uedge import bbb
        for key, value in self.orig.items():
            bbb.__setattr__(key, value)

    def exmain_isaborted(self):
        ''' Checks if abort is requested '''
        from uedge import bbb
        bbb.exmain()
        # Abort flag set, abort case
        if bbb.exmain_aborted == 1: 
            # Reset flag
            bbb.exmain_aborted == 0 
            # Restore parameters modified by script
            try:
                self.restorevalues(self)
            except:
                self.restorevalues()
            return True

    def message(self, string, separator='-', pad='', seppad = '', 
        nseparator=1):
        ''' Prints formatted self.message to stdout '''
        # TODO: add formatting for len>75 strings
        msg = pad.strip() + ' ' + string.strip() + ' ' + pad.strip()
        for i in range(nseparator):
            print(seppad + separator*(len(msg)-2*len(seppad)) + seppad)
        print(msg)
        print(seppad + separator*(len(msg)-2*len(seppad)) + seppad)




    def converge(self, dtreal=2e-9, ii1max=5000, ii2max=5, itermx=7, ftol=1e-5,
        dt_kill=1e-14, t_stop=100, dt_max=100, ftol_min = 1e-9, incpset=7,
        n_stor=0, storedist='lin', numrevjmax=2, numfwdjmax=1, numtotjmax=0, 
        tstor=(1e-3, 4e-2), ismfnkauto=True, dtmfnk3=5e-4, mult_dt=3.4, 
        reset=True, initjac=True, rdtphidtr=1e20, deldt_min=0.04, rlx=0.9,

        tsnapshot=None, savedir='../solutions', ii2increase=0, savefname=None,
        message=None):

        ''' Converges the case by increasing dt 
        dtreal : float [1e-9]
            Original time-step size
        ii1max : int [500]
            Outer loop iterations, i.e. time-step changes
        ii2max : int [5]
            Inner loop iterations, i.e. time-steps per time-step change
        dt_kill : float [1e-14]
            Time-step limit for aborting simulation
        itermx : int [7]
            Maximum iterations per time-step used internally in routine
        ftol : float [1e-5]
            Internal fnrm tolerance for time-steps
        incpset : int [7]
            
        savedir : str ['../solutions']

        numtotjmax : int [None]

        ii2increase : float [1.5]
            
        ftol_min : float [1e-9]
            Value of fnrm where time-advance will stop
        t_stop : float [100.]
            Maximum total accumulated plasma-time before stopping if 
            fnorm has not decreased below ftol_min
        dt_max : float [100.]
            Maximum allowable time-step size
        numrevjmax : int [2]
            Number of time-step reducitons before Jacobian is recomputed
        numfwdjmax : int [2]
            Number of time-step increases before Jacobian is recomputed
        n_stor : int [0]
            Number of time-slices to be stored in interval tstor
        tstor : tuple of floats [(1e-3, 4e-2)]
            Time-interval in which time-slices are stored (lower, upper)
        storedist : str ['lin']
            Distribution of time-slices in tstor. Options are 'lin' and 
            'log' for linear and logarithmic distributions, respectively
        reset : bool [True]
            Switch whether to reset the total time etc of the case
        initjac : bool [False]
            Switch to re-evaluate Jacobian on first convegnec time-step
            or not
        ismfnkauto : bool [True]
            If True, sets mfnksol=3 for time-steps smaller that dtmfnk3,
            mfnksol=-3 for larger time-step sizes
        dtmfnk3 : float [5e-4]
            Time-step size for which ismfnkauto controls mfnksol if
            ismfnkauto is True
        mult_dt : float [3.4]
            Time-step size increase factor after successful inner loop
        rdtphidtr : float [1e20]
            Ratio of potential-equation time-step to plasma equaiton time-step
            size: dtphi/dtreal
        deldt_min : float [0.04]
            Minimum relative change allowed for model_dt>0
        rlx : float [0.9]
            Maximum change in variable at each internal linear iteration
        tsnapshot : list [None]
            If None, uses linear/logarithmic interpolation according to 
            storedist in the interval tstor. Snapshot times can be defined in 
            a list and supplied. Then, the tnsaphsot list defines the 
            time-slices
         
        '''
        from numpy import linspace, logspace, log10, append
        from copy import deepcopy
        from uedge import bbb
        from time import time
        from os.path import exists

        self.tstart = time()
        # TODO: add kwarg savename

        # Check if requested save-directory exists: if not, write to cwd
        if not exists(savedir):
            print('Requested save-path {} not found, writing to cwd!'.format(\
                savedir))
            savedir = '.'        

        if savefname is None:
            autoname = True
            try:
                savefname = self.casename
            except:
                savefname = 'dtrun'
        else:
            autoname = False
        
        savefname = '{}/{}'.format(savedir, savefname)

        if autoname is True:
            if exists('{}_last_ii2.hdf5'.format(savefname)):
                i = 2
                while exists('{}_{}_last_ii2.hdf5'.format(savefname, i)):
                    i += 1
                savefname = '{}_{}'.format(savefname, i)
        savefname = '{}_last_ii2.hdf5'.format(savefname)
             
        self.classvars = {}
        for  var in ['slice_ni', 'slice_ng', 
            'slice_up', 'slice_te', 'slice_ti', 'slice_tg', 'slice_phi', 
            'slice_dt_tot', 'time', 'fnorm', 'nfe', 'dt_tot', 'dtreal', 'ii1', 
            'ii2', 'ii1fail', 'ii2fail', 'dtrealfail', 'itrouble', 'troubleeq',
            'troubleindex', 'ylfail', 'internaleq', 'internalspecies', 
            'yldotsfscalfail']:
            self.classvars[var] = []

        self.orig = {}
        for var in ['itermx', 'dtreal', 'icntnunk', 'ftol', 'mfnksol', 'rlx',
            'deldt', 'isdtsfscal', 'incpset']:
            self.orig[var] = deepcopy(getattr(bbb, var))
        if numtotjmax == 0:
            numtotjmax = numrevjmax + numfwdjmax
        self.classsetup = {}
        for var in ['itermx', 'incpset', 'ii1max', 'ii2max', 'numrevjmax',
            'numfwdjmax', 'numtotjmax', 'rdtphidtr', 'deldt_min', 'rlx']:
            self.classsetup[var] = locals()[var]
# TODO: Add variable to control reduciton factor?
# TODO: Should dtreal = min(x, t_stop) actually be t_stop or dt_max?

        def scale_timestep(scaling):
            ''' Increases/decreases time-step ''' 
            bbb.dtreal *= scaling

        def issuccess(self, t_stop, ftol_min, savefname):
            ''' Checks if case is converged '''
            from datetime import timedelta
            from time import time
            if (bbb.iterm == 1):
                bbb.ylodt = bbb.yl
                bbb.dt_tot += bbb.dtreal
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                self.fnrm_old = sum((bbb.yldot[:bbb.neq-1]*\
                    bbb.sfscal[:bbb.neq-1])**2)**0.5
                self.savesuccess(savefname,
                    #bbb.label[0].strip(\).decode('UTF-8'), 
                    self.fnrm_old)
                if (bbb.dt_tot>=t_stop  or  self.fnrm_old<ftol_min):
                    print('')
                    self.message('SUCCESS: ' + 'fnrm < bbb.ftol'\
                        *(self.fnrm_old<ftol_min) + \
                        'dt_tot >= t_stop'*(bbb.dt_tot >= t_stop), pad='**', 
                        separator='*')
                    print('Total runtime: {}'.format(timedelta(
                        seconds=round(time()-self.tstart))))
                    try:
                        self.restorevalues(self)
                    except:
                        self.restorevalues()
                    return True
    
        def isfail(dt_kill):
            ''' Checks whether to abandon case '''
            if (bbb.dtreal < dt_kill):
                self.message('FAILURE: time-step < dt_kill', pad='**', 
                separator='*')
                try:
                    self.restorevalues(self)
                except:
                    self.restorevalues()
                return True

        def setmfnksol(ismfnkauto, dtmfnk3):
            ''' Sets mfnksol according to setup '''
            if ismfnkauto is True:
                bbb.mfnksol = 3*(-1)**(bbb.dtreal > dtmfnk3)

        def calc_fnrm():
            ''' Calculates the initial fnrm '''
            from uedge import bbb
            bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
            return sum((bbb.yldot[:bbb.neq-1]*bbb.sfscal[:bbb.neq-1])**2)**0.5
 
        ''' TIME-SLICING SETUP '''
        if tsnapshot is None:
            if storedist == 'lin':
                # Linearly spaced time slices for writing 
                dt_stor = linspace(tstor[0], tstor[1], n_stor)
            elif storedist == 'log':
                # Logarithmically spaced time-slices
                dt_stor = logspace(log10(tstor[0]), log10(tstor[1]), n_stor)
        else:
            dt_stor = tsnapshot
        # Add end-point to avoid tripping on empty arrays
        dt_stor = append(dt_stor, 1e20)

        if reset is True:
            bbb.dt_tot = 0               

        ''' TIME-STEP INITIALIZATION '''
        bbb.rlx = rlx
        bbb.dtreal = dtreal
        bbb.ftol = ftol
        ii1 = 0
        ii2 = 0
        if (bbb.iterm == 1) and (bbb.ijactot > 0):
            self.message('Initial successful time-step exists', separator='')
        else:
            self.message('Need to take initial step with Jacobian; ' + \
                'trying to do here', seppad='*')
            # Ensure time-step is taken
            bbb.icntnunk = 0
            if message is not None:
                print(message)
            # Take timestep and see if abort requested
            if self.exmain_isaborted():
                return
            # Increase time
            # Verify time-step was successful
            if (bbb.iterm != 1):
                try:
                    self.restorevalues(self)
                except:
                    self.restorevalues()
                self.message('Error: converge an initial time-step first; '
                    'then re-execute command', seppad='*')
                return
        bbb.incpset = incpset
        bbb.itermx = itermx
        deldt_0 = deepcopy(bbb.deldt)
        isdtsf_sav = deepcopy(bbb.isdtsfscal)
# TODO: Replace with some more useful information?
#        if (bbb.ipt==1 and bbb.isteon==1): # set ipt to te(nx,iysptrx+1) 
#           #if no user value
#           ipt = bbb.idxte[nx-1,com.iysptrx] #note: ipt is local, 
#               # bbb.ipt global
        bbb.dtphi = rdtphidtr*bbb.dtreal
        svrpkg=bbb.svrpkg.tostring().strip()
        bbb.ylodt = bbb.yl
        self.fnrm_old = calc_fnrm()
        if initjac is True: 
            self.fnrm_old = 1e20
        else:
            bbb.newgeo=0
        # Intialize counters
        irev = -1
        numfwd = 0
        numrev = 0
        numrfcum = 0
        # Compensate for first time-step before entering loop
        scale_timestep(1/(3*(irev == 0) + mult_dt*(irev != 0)))
        ''' OUTER LOOP - MODIFY TIME-STEP SIZE'''
        # TODO: Add logic to always go back to last successful ii2 to 
        # precondition the Jacobian, to avoid downwards cascades?
        # NOTE: experomental functionality
        successivesuccesses = 0
        for ii1 in range(ii1max):
            setmfnksol(ismfnkauto, dtmfnk3)
            # adjust the time-step
            # dtmult=3 only used after a dt reduc. success. completes loop ii2
            # for fixed dt either increase or decrease dtreal; depends 
            # on mult_dt
            scale_timestep(3*(irev == 0) + mult_dt*(irev != 0))
            bbb.dtreal = min([bbb.dtreal, dt_max]) 
            bbb.dtphi = rdtphidtr*bbb.dtreal
            bbb.deldt =  min([bbb.deldt, deldt_0, deldt_min])
            self.message('Number of time-step changes = ''{} New time-step: {:.2E}\n'\
                .format((ii1+1), bbb.dtreal), pad='***', nseparator=1)

            # Enter for every loop except first, unless intijac == True
            if ii1 > -int(initjac): 
                # Main time-stepping switch: controls increase/decrease in 
                # dtreal and Jacobian preconditioning
                if (irev == 1):      # decrease in bbb.dtreal
                    if (numrev < numrevjmax and \
                        numrfcum < numtotjmax): #dont recom bbb.jac
                        bbb.icntnunk = 1	
                        numrfcum += 1
                    else:               # force bbb.jac calc, reset numrev
                        bbb.icntnunk = 0
                        numrev = -1		# yields api.zero in next statement
                        numrfcum = 0
                    numrev += 1
                    numfwd = 0
                else:  # increase in bbb.dtreal
                    if (numfwd < numfwdjmax and \
                        numrfcum < numtotjmax): 	#dont recomp bbb.jac
                        bbb.icntnunk = 1
                        numrfcum += 1
                    else:
                        bbb.icntnunk = 0 # recompute jacobian for increase dt
                        numfwd = -1
                        numrfcum = 0
                    numfwd += 1
                    numrev = 0			#bbb.restart counter for dt reversals
                bbb.isdtsfscal = isdtsf_sav
                # Dynamically decrease ftol as the initial ftol decreases
                bbb.ftol = max(min(ftol, 0.01*self.fnrm_old),ftol_min)
                # Take timestep and see if abort requested
                if message is not None:
                    print(message)
                if self.exmain_isaborted():
                    return
                if bbb.iterm == 1:
                    self.classvars['ii1'].append(ii1)
                    self.classvars['ii2'].append(ii2)
                if issuccess(self, t_stop, ftol_min, savefname):
                    return
            bbb.icntnunk = 2
            bbb.isdtsfscal = 0
            # NOTE: experomental functionality
            bbb.ii2max = ii2max + round(ii2increase*successivesuccesses)
            # Take ii2max time-steps at current time-step size while 
            # time-steps converge: if not, drop through
            for ii2 in range(bbb.ii2max): 
                if (bbb.iterm == 1):
                    bbb.ftol = max(min(ftol, 0.01*self.fnrm_old),ftol_min)
                    # Take timestep and see if abort requested
                    self.message("Inner iteration #{}".format(ii2+1), nseparator=0, 
                        separator='')
                    if message is not None:
                        print(message)
                    if self.exmain_isaborted():
                        return
                    if bbb.iterm == 1:
                        self.classvars['ii1'].append(ii1)
                        self.classvars['ii2'].append(ii2)
                    if issuccess(self, t_stop, ftol_min, savefname):
                        return
                    self.message("Total time = {:.4E}; Timestep = {:.4E}".format(\
                        bbb.dt_tot-bbb.dtreal,bbb.dtreal), nseparator=0, 
                        separator='')
# TODO: replace with more useful information
#                       print("variable index ipt = ",ipt, " bbb.yl[ipt] = ", 
#                            bbb.yl[ipt])
                    # Store variable if threshold has been passed
                    if (bbb.dt_tot >= dt_stor[0]):
                        # Remove storing time-points smaller than current 
                        # simulation time
                        while bbb.dt_tot >= dt_stor[0]:
                            dt_stor = dt_stor[1:]
                        self.store_timeslice()
            irev -= 1
            # Output and store troublemaker info
            # NOTE: experomental functionality
            successivesuccesses += 1
            if (bbb.iterm != 1):	
                # NOTE: experomental functionality
                successivesuccesses = 0
                self.itroub()
                ''' ISFAIL '''
                if isfail(dt_kill):
                    self.save_intermediate(savefname)
#                        bbb.label[0].strip().decode('UTF-8'))
                    break
                irev = 1
                self.message('Converg. fails for bbb.dtreal; reduce time-step by '+\
                    '3, try again', pad = '***', nseparator=0)
                scale_timestep(1/(3*mult_dt))
                bbb.dtphi = rdtphidtr*bbb.dtreal
                bbb.deldt *=  1/(3*mult_dt) 
                setmfnksol(ismfnkauto, dtmfnk3)
#                bbb.iterm = 1
        # ii1max iterations completed without convergence
        # Indicate solver failure regardless of whether last solve was completer
        # or not.
        bbb.iterm = -1 # Ensure subsequent repetitions work as intended
        # Return False to signal fail
        return False


    def continuation_solve(self, 
        var, 
        target=None, 
        savedir=None,
        index=None, 
        dt=0.02,
        dtreal=1e-5, 
        ftol=1e-5,
        iicont_max=7,
        iicont_fail_max=3,
        cutoff=1e6,
        itermx=5,
        incpset=8,
        initres=1000,
        deltafac=2,
        dtdeltafac=2,
        staticiter_max=3,
        dtlim=1e4,
        ii1max=150,
        saveres=7,
        newgeo=False,
        commands = [],
        **kwargs):
        """ Solves var up to target using continuation method
        var - string for variable, or nested dict of variables and targets
                If nested dict, top level is variable name and subdict contains
                keys 'target' and optionally 'index'. These values override the 
                kwargs of the solver. If index is not set, defaults to None
        target - target value for variable
        index - tuple containing index of var if var is array.
                In case of extended array slicing, e.g. setting ranges of
                var, index must be defined as tuples of indices and/or 
                slice objects. E.g. to set var[1:4,:,0], one would use
                index = (slice(1,4), slice(None), 0). If ranges of var is
                set, ensure target has the matching, appropriate dimensions
                or is a float
        dt - time-step used when solving using continuation
        staticiter_max - maximum consequtive static iterations before callong 
                time-depenent run
        dtreal - time-step used when resorting to time-dependent solve
        ftol - ftol to be attained
        iicont_max - iterations in continuation loop
        iicont_fail_max - max fails in continuation loop before going time-dependent
        cutoff - cutoff value for required interations at current delta to attain target
        itermx - max number of iterations to take within time-step
        incpset - max iterations before re-calculating jacobian
        initres - inital resolution, i.e. delta = target-var/initres
        deltafac - factor by which delta is increased after iicont_max increases
        dtdeltafac - factor by which delta is increased by for time-dependent run
        dtlim - Cutoff on (target-var)/delta for which a time-dependent run is triggered
        ii1max - maximum number of dt iterations to take
        saveres - how many successful iterations are allowed between saves are dumped
        newgeo - UEDGE re-evaluates the grid if newgeo is True
        kwargs passed to rundt
        """
        from uedge import bbb
        from time import time
        from os import makedirs
        from copy import deepcopy
        from shutil import rmtree
        from os.path import exists
        from numpy import array
        from Forthon import package, packageobject
        # TODO: icntnunk=0 on fail only?

        self.tstart = time()
        # ==== HELPERS ====
        def changevar(commands=[]):
            """ Changes var by delta """
            setvar(min(self.lastsuccess + self.delta, 1))
            for command in self.contsolvecommands:
                exec(command)
           
        def setvar(value):
            """ Sets var to value """
            for key, subdict in self.var.items():
                newvar = subdict['origvar'] + value * subdict['deltavar']
                if subdict['index'] is None:
                    setattr(subdict['pkgobj'], key, newvar)
                else:
                    getattr(subdict['pkgobj'], key)[subdict['index']] = newvar

        def getvar(key, subdict):
            """ Returns current value of var """
            # Variable is array: only set requested index
            if 'index' not in subdict:
                return getattr(subdict['pkgobj'], key)
            else:
                try:
                    return getattr(subdict['pkgobj'], key)[subdict['index']]
                except:
                    return getattr(subdict['pkgobj'], key)

        def issuccess():
            """ Checks whether case is converged at target value of var """
            from time import time
            from datetime import timedelta
            self.lastsuccess += min(self.delta, 1-self.lastsuccess)
            if self.lastsuccess >= 1:
                print('\n===== TARGET VALUE ACHIEVED: REDUCE FNORM ====')
                staticiter()
                bbb.dtreal = 1e20

                print('++++++++++++++++++++++++++++++++++++++++')
                print('+++++ CONTINUATION SOLVE SUCCEEDED +++++')
                print('++++++++++++++++++++++++++++++++++++++++')
                print('Total runtime: {}'.format(timedelta(
                        seconds=round(time()-self.tstart))))
                self.savesuccess('{}/SUCCESS.hdf5'.format(self.savedir))
                self.restorevalues()
                return True
            else:
                if self.isave >= self.saveres:
                    self.isave = 0
                    self.savesuccess(self.savefname.format('{:.3f}'.format(100*self.lastsuccess).replace('.','p')))
                else:
                    self.isave += 1

        def printinfo(printdelta=True):
            """ Print current solve status """
            # TODO: add printing of indices?
            ljust = 42
            print('    -Variable(s) being solved:')
            for key, _ in self.var.items():
                print(' '.ljust(ljust-3), '-',key)
            print('{}{:.3f}%'.format('    -Progress'.ljust(ljust),
                        self.lastsuccess*100
            ))
            if printdelta is not False:
                print('{}{:.3f}%'.format('    -Advancing by'.ljust(ljust),
                        self.delta*100
                ))
                print('{}{}'.format('    -Steps to target at current delta:'.ljust(ljust), 
                    (int((1 - self.lastsuccess)/self.delta))))
            print('')

        def isdeltatoosmall():
            """ Returns True if delta is below cutoff """
            if self.delta < 1/self.cutoff:
                return True
            return False

        def dostaticiter():
            """ Returns True if delta is above dt-run threshold """
            if self.delta > 1/self.dtlim:
                return True
            return False
                    
        def staticiter():
            """ Tries to reduce initial fnrm by iterating at constant var """
            from uedge import bbb
            from copy import deepcopy
            # TODO: always dump save on successful steady-state solution?
            printinfo(False)
            setvar(self.lastsuccess)
            dtreal_orig = deepcopy(bbb.dtreal)
            ftol_orig = deepcopy(bbb.ftol)
            bbb.dtreal = dtreal_orig/100
            while bbb.dtreal < 1e5:
                ftol_old = deepcopy(bbb.ftol)
                # Dynamically decreasing fnorm
                bbb.pandf1 (-1, -1, 0, bbb.neq, 1., bbb.yl, bbb.yldot)
                fnorm_old = (sum((bbb.yldot[:bbb.neq]*\
                    bbb.sfscal[:bbb.neq])**2))**0.5
                if bbb.dtreal > 1:
                    bbb.dtreal = 1e20
                print('\n===== STATIC ITERATION FOR DTREAL={:.2e} AT {:.1f}% ====='.format(bbb.dtreal, self.lastsuccess*100))
                bbb.ftol = max(min(bbb.ftol, 0.01*fnorm_old), 1e-9)
                # Take a converging step
                if self.exmain_isaborted():
                    setvar(self.lastsuccess)
                    bbb.dtreal = dtreal_orig
                    bbb.ftol = ftol_orig
                    return None
                # Temporary fix to avoid Segfault when using icntnunk=1!
                bbb.ijactot = 2
                bbb.icntnunk = 1
                bbb.ftol = ftol_old
                bbb.dtreal = bbb.dtreal*5
                # Assert steady-state
                if bbb.iterm != 1:
                    bbb.dtreal = dtreal_orig
                    setvar(self.lastsuccess)
                    bbb.ftol = ftol_orig
                    print('\n===== STATIC FNRM REDUCTION FAILED =====\n')
                    return False
            self.savesuccess(self.savefname.format('{:.3f}_staticiter'.format(\
                    100*self.lastsuccess).replace('.','p')
                ))
            print('===== CONVERGED AT STEADY STATE: RETURNING TO MAIN LOOP =====')
            bbb.dtreal = dtreal_orig
            bbb.ftol = ftol_orig
            return True

        def dtsolve(dtdeltafac):
            """ Perform a time-dependent solve """
            from uedge import bbb
            printinfo(False)
            bbb.icntnunk = 0
            # Make backup of original values before entering dt run
            # TODO: compress the storing/setting/resetting of vars
            tstart_cont = deepcopy(self.tstart)
            incpset_cont = deepcopy(bbb.incpset)
            itermx_cont = deepcopy(bbb.itermx)
            dtreal_cont = deepcopy(bbb.dtreal)
            ftol_cont = deepcopy(bbb.ftol)
            orig_cont = deepcopy(self.orig)
            classvars_cont = deepcopy(self.classvars)
            classsetup_cont = deepcopy(self.classsetup)
            bbb.incpset = 5
            bbb.itermx = 30
            bbb.dtreal = 1e-5
            bbb.ftol = 1e-8
            abort = False
            # Ensure a first time-step can be taken
            dtdelta = self.lastsuccess + dtdeltafac/100
            while bbb.iterm != 1:
                dtdelta = min(self.lastsuccess + dtdeltafac/100, 1)
                setvar(dtdelta)
                if self.exmain_isaborted():
                    setvar(self.lastsuccess)
                    return False
                # If the iteration does not succeed, reduce the time-step
                if bbb.iterm != 1:
                    # Advancing less than .5% dt is probably an indication of a failure
                    if dtdeltafac < 0.5:
                        abort = True
                        break
                    dtdeltafac /=2
            if abort is False:
                self.converge(dtreal=dtreal, savedir='.', 
                    savefname=self.savefname.format('{:.3f}_dtrun'.format(\
                    dtdelta*100).replace('.','p')).replace('.hdf5',''), 
                    message='Solving for delta={:.3f}%'.format(dtdelta*100),
                    ii1max=ii1max, **kwargs)
                if bbb.iterm == 1:
                    self.lastsuccess = dtdelta
                # Ensure original values are still being kept track of
                self.orig = orig_cont
                bbb.incpset = incpset_cont 
                bbb.itermx = itermx_cont 
                bbb.dtreal = dtreal_cont 
                bbb.ftol = ftol_cont 
                orig_cont = deepcopy(self.orig)
                self.classvars = classvars_cont
                self.classsetup = classsetup_cont
                self.tstart = tstart_cont
                bbb.dtreal = dt
            if bbb.iterm != 1:
                print('=====================================')
                print('===== CONTINUATION SOLVE FAILED =====')
                print('=====================================')
                setvar(self.lastsuccess)
                print('Progress upon abortion: {:.3f}%'.format(
                    self.lastsuccess*100)
                )
                return False
            else:
                return True

        if isinstance(var, str): 
            if target is None:
                raise KeyError("Keyword argument 'target' must be set"
                    " when evaluating a single variable!"
                )
        elif isinstance(var, dict):
            if savedir is None: 
                raise KeyError("Must define save directory name when "
                    "evalauting multiple variables!"
                )
        self.savedir = savedir
        if self.savedir is None:
            self.savedir = bbb.label[0].decode('UTF-8').strip()
            autosavedir = True
        else:
            autosavedir = False
        if len(self.savedir) == 0:
            self.savedir = '{}_solve'.format(var)
            autosavedir = True
        if autosavedir:
            if exists(self.savedir):
                i=2
                while exists(self.savedir+'_{}'.format(i)):
                    i += 1
                print('Run-folder {} exists: saving to {}_{}'.format(\
                    self.savedir, self.savedir, i)
                )
                self.savedir = '{}_{}'.format(self.savedir, i)
            else:
                print('Saving intermediate saves to {}'.format(self.savedir))
        else:
            if exists(self.savedir):
                rmtree(self.savedir)
        # Finalize savename w/ placeholder
        self.savefname = '{}/progress{{}}.hdf5'.format(self.savedir)
        makedirs(self.savedir)
        # Additional variables to be stored in save files
        self.classvars = {}
        for dictkeys in ['delta', 'progress', 'delta_fail', 'progress_fail']:
            self.classvars[dictkeys] = []
        # Find package variable
        # TODO: Use dict to set one or multiple variables instead of string?

        # Set up variables that need to be accessed by helpers
        self.dtlim = dtlim
        self.cutoff = cutoff
        # Set up control flags
        start = True
        nstaticiter = 0
        iicont = 0
        iicont_fail = 0
        dtiter = 0
        self.isave = saveres
        self.saveres = saveres
        self.delta = 1/initres
        self.contsolvecommands = commands
        self.lastsuccess = 0

        # Set up dictionary with variables and targets
        self.var = {}
        # Construct the variable dictionary based on the user input
        if isinstance(var, str):
            self.var = {    var: {
                            'index': index,
                            'target': target
                    }}
        else:
            self.var = var
            for key, subdict in self.var.items():
                if 'target' not in subdict.keys():
                    raise KeyError("Target not set for variable {}!".format(key))
                elif 'index' not in subdict.keys():
                    self.var[key]['index'] = None
        # Complement the user input with gradients for each variable
        for key, subdict in self.var.items():
            if isinstance(subdict['target'], list):
                subdict['target'] = array(subdict['target'])
            for pkg in package():
                pkgobj = packageobject(pkg)
                if key in pkgobj.varlist():
                    subdict['pkg'] = pkg
                    subdict['pkgobj'] = pkgobj
            if subdict['index'] is not None:
                try:
                    getattr(subdict['pkgobj'], key)[subdict['index']]
                except:
                    raise ValueError('{} does not accommodate requested indices!'.format(key))
            subdict['origvar'] = deepcopy(getvar(key, subdict))
            subdict['deltavar'] = deepcopy(subdict['target'] - subdict['origvar'])
            try:
                dist = abs(subdict['deltavar']).max()
            except:
                dist = abs(subdict['deltavar'])
            if dist < 1e-10:
                raise ValueError('Target equals current value for {}, aborting!'.format(key))
                
        self.classsetup = {}
        for key, subdict in self.var.items():
            self.setupvars[key] = getattr(subdict['pkgobj'], key)
            self.classsetup['initial_{}'.format(key)] = getvar(key, subdict)
            self.classsetup['delta_{}'.format(key)] = subdict['deltavar']
            self.classsetup['target_{}'.format(key)] = subdict['target']
            if isinstance(subdict['index'], (tuple, slice)):
                self.classsetup['index_{}'.format(key)] = str(subdict['index']).encode("ascii", "ignore")
            elif subdict['index'] is not None:
                self.classsetup['index_{}'.format(key)] = subdict['index']
            else:
                self.classsetup['index_{}'.format(key)] = False

        # Record original solver settings            
        self.orig = {}
        for ovar in ['itermx', 'dtreal', 'icntnunk', 'ftol', 'incpset',
            'ismmaxuc', 'mmaxu', 'newgeo' 
        ]:
            self.orig[ovar] = deepcopy(getattr(bbb, ovar))
        # Take the initial time-step
        changevar()
        # Set remaining solver settings
        bbb.dtreal = dt
        bbb.ftol = ftol
        bbb.incpset = incpset
        bbb.itermx = itermx
        # TODO: Resolve how to run continuation solver w/ mmaxu
        # TODO: Add coding for autodetecting when avg nfe approx mmaxu
        #       and upate jacobian then

#        bbb.ismmaxuc = 0
#        bbb.mmaxu = 70
        if (bbb.iterm == 1) and (bbb.ijactot > 0):
            self.message('Initial successful time-step exists', separator='')
        else:
            self.message('Need to take initial step with Jacobian; ' + \
                'trying to do here', seppad='*')
            printinfo()
            # Ensure preconditioner is calculated
            bbb.icntnunk = 0
            # Take timestep and see if abort requested
            if self.exmain_isaborted():
                setvar(self.lastsuccess)
                return
            # Increase time
            # Verify time-step was successful
            if (bbb.iterm != 1):
                self.restorevalues()
                setvar(self.lastsuccess)
                self.message('Error: converge an initial time-step first; then ' + \
                    're-execute command', seppad='*')
                return
        # Start outer loop

        if newgeo is False:
            bbb.newgeo = False

        changevar()
        while True:
            # Ensure we don't exceed the target value
            # Re-eval preconditioner: omit first step
            if start is False:
                bbb.icntnunk = 0
            else: 
                bbb.icntnunk = 1
                start = False
            # Reset counters
            iicont = 0
            iicont_fail = 0
            # Start inner loop
            while iicont < iicont_max:
                print('===== MAIN LOOP {}/{} ====='.format(iicont+1, iicont_max))
                printinfo()
                # TODO: how to get intial fnrm rather than last fnorm?
                fnorm_old = (sum((bbb.yldot[:bbb.neq]*\
                    bbb.sfscal[:bbb.neq])**2))**0.5
                bbb.ftol = min(max(fnorm_old/100, 1e-6), ftol)
                # Take step and check whether abort is requested
                if self.exmain_isaborted():
                    setvar(self.lastsuccess)
                    return
                bbb.ftol = ftol
                bbb.icntnunk = 1
                # Failure to take time-step
                if bbb.iterm != 1:
                    # Check whether an excorbant number of iterations are necessary to converge: fail if yes
                    self.classvars['delta_fail'].append(self.delta)
                    self.classvars['progress_fail'].append(self.lastsuccess)
                    if isdeltatoosmall():
                        print('=====================================')
                        print('===== CONTINUATION SOLVE FAILED =====')
                        print('=====================================')
                        print('Last successful step for {}: {:.4e}'.format(\
                            str(list(self.var.keys())).replace('[','').replace(']',''), 
                            self.lastsuccess)
                        )
                        return
                    # Start trying to reset convergence
                    else:
                        # Reset success counter
                        iicont = 0
                        iicont_fail += 1
                        # Delta has been reduced too many times: do some moves to knock the case loose
                        if iicont_fail >= iicont_fail_max:
                            print('===== RECOVERY LOOP ENTERED =====')
                            # Iterate the case at the last successful 
                            bbb.icntnunk = 0
                            if dostaticiter():
                                staticstatus = staticiter()
                                # Flag for abortion
                                if staticstatus is None:
                                    return 
                                elif staticstatus is False:
                                    dtcall=True
                                else:
                                    dtcall=False
                                nstaticiter += 1
                            # Enter loop for time-dependent simulations
                            if (bbb.iterm != 1) or (dtcall is True) or (nstaticiter==staticiter_max):
                                bbb.iterm = 0
                                print('===== ENTER TIME-DEPENDENT SOLVE =====')
                                if dtsolve(dtdeltafac) is False:
                                    return
                                else:
                                    # Force-break inner loop
                                    iicont = 1e5
                                    dtcall = False
                                    bbb.icntnunk = 0
                                    if issuccess():
                                        return
                                    changevar()
                            # Initial fnrm successfully reduced, continue increasing var
                            else:
                                iicont_fail = 0
                                changevar()
                        # Try to attain convergence again by decreasing change in var
                        else:
                            print('\n===== FAIL {}/{} ====='.format(iicont_fail, iicont_fail_max))
                            print('===== DECREASE DELTA AND TRY AGAIN =====')
                            # Reset to last successful step
                            setvar(self.lastsuccess)
                            # Decrease change in variable
                            self.delta /= deltafac
                            # Set to new, decreased change
                            changevar()
                            bbb.icntnunk = 0
                # Success: update counter and last successful value
                else:
                    self.classvars['delta'].append(self.delta)
                    self.classvars['progress'].append(self.lastsuccess)
                    if issuccess():
                        return
                    # Check whether to increment time-step
                    if iicont == iicont_max-1:
                        print('\n===== INNER LOOP COMPLETED: ADVANCING DELTA =====')
                        nstaticiter = 0
                        self.delta *= 1.1*deltafac
                    else:
                        print('\n===== SUCCESS: ADVANCING VARIABLE =====')
                    iicont += 1
                    changevar()
        print('EXITED LOOP: YOU SHOULD NOT BE HERE...')
        return

def rundt(**kwargs):
    runcase=UeRun()
    runcase.converge(**kwargs)
