
def write_subpandf1():
    lines = write_ompsubroutine('convsr_vo1', "(xc, yc, ylcopy)", True)
    lines = write_ompsubroutine('convsr_vo2', "(xc, yc, ylcopy)", True)
    lines = write_ompsubroutine('convsr_aux1', "(xc, yc)", True)
    lines = write_ompsubroutine('convsr_aux2', "(xc, yc)", True)
    lines = write_ompsubroutine('calc_plasma_diffusivities', "", True)
    lines = write_ompsubroutine('initialize_driftterms', "", True)
    lines = write_ompsubroutine('calc_driftterms1', "", True)
    lines = write_ompsubroutine('calc_driftterms2', "", True)
    lines = write_ompsubroutine('calc_currents', "", True)
    lines = write_ompsubroutine('calc_fqp', "", True)
    lines = write_ompsubroutine('calc_friction', "(xc)", True)
    lines = write_ompsubroutine('calc_elec_velocities', "", True)
    lines = write_ompsubroutine('calc_volumetric_sources', "(xc,yc)", True)
    for line in lines:
        print(line)
    

def write_ompsubroutine(subroutine, arguments, bounds=False, subcalls=[]):
    from uetools import Case
    from textwrap import wrap

    outlines = []
    getpkg = Case().search.get_group
    pkgs = {}
    defined_vars = get_pandf_vars()
    varlist = []
    varlist = varlist + defined_vars[subroutine]
    for subcall in subcalls:
        varlist = varlist + defined_vars[subcall]
    for var in varlist:
        pkg = getpkg(var)
        if pkg not in pkgs:
            pkgs[pkg] = []
        pkgs[pkg].append(var)
    tmpvarlist = [f"{x}_tmp" for x in varlist]

    bbbv = parse_bbbv()
    vardef = []
    for var in varlist:
        vardef.append("{}_tmp({})".format(var, ",".join(bbbv[var])))



    outlines.append(f"  SUBROUTINE OMP{subroutine}(neq, yl, yldot)")
    outlines.append( "    USE Dim, ONLY: nx, ny, ngsp, nisp, nxpt")
    outlines.append( "    USE OMPPandf1Settings, ONLY:OMPPandf1loopNchunk")
    outlines.append( "    USE OmpCopybbb")


    for pkg, var in pkgs.items():
        lines =  [x.replace(" ",", ") for x in wrap( 
                    " ".join(var), width=70, break_long_words=False)]
        outlines.append("    USE {}, ONLY: {}".format(pkg, ", &\n    &    ".join(lines)))
    outlines.append("    IMPLICIT NONE")
    outlines.append("    INTEGER, INTENT(IN):: neq")
    outlines.append("    REAL, INTENT(IN):: yl(*)")
    outlines.append("    REAL, INTENT(OUT):: yldot(*)")
    outlines.append("    INTEGER:: chunks(1:neq,3), Nchunks, ichunk, xc, yc")
    outlines.append("    REAL:: yldotcopy(1:neq), ylcopy(1:neq+2)")
    outlines.append("! Define local variables")
    lines =  [x.replace(" ",", ") for x in wrap(
            " ".join(vardef), width=70, break_long_words=False)]
    outlines.append("    real:: {}".format(", &\n    &      ".join(lines)))

    outlines.append("\n    ! Initialize arrays to zero")        
    lines =  [x.replace(" ","; ") for x in wrap( 
                " ".join([f'{x}=0.' for x in tmpvarlist]), 
                width=70, break_long_words=False)]
    outlines.append(4*" " + f"\n    ".join(lines))
    outlines.append("\n    ylcopy(1:neq+1)=yl(1:neq+1); yldotcopy=0")
    
    if bounds:
        outlines.append("\n    call chunk3d(0,nx+1,0,ny+1,0,0,chunks,Nchunks)")
    else:
        outlines.append("\n    call chunk3d(1,nx,1,ny,0,0,chunks,Nchunks)")


    outlines.append("\n    !$OMP    PARALLEL DO &")
    outlines.append("    !$OMP &      default(shared) &")
    outlines.append("    !$OMP &      schedule(dynamic,OMPPandf1LoopNchunk) &")
    outlines.append("    !$OMP &      private(ichunk,xc,yc) &")
    outlines.append("    !$OMP &      firstprivate(ylcopy, yldotcopy) &")
    lines =  [x.replace(" ",", ") for x in wrap( 
                " ".join(tmpvarlist), width=70, break_long_words=False)]
    outlines.append(4*" " + "!$OMP &      REDUCTION(+:{})".format(\
            ', &\n    !$OMP &         '.join(lines)))


    outlines.append(4*" " + "DO ichunk = 1, Nchunks")
    outlines.append(8*" "+"xc = chunks(ichunk,1)")
    outlines.append(8*" "+"yc = chunks(ichunk,2)")
    outlines.append(8*" "+"call initialize_ranges(xc, yc, 0, 0, 0)")
    outlines.append(8*" " + "call {}{}".format(subroutine, arguments))


    outlines.append("\n        ! Update locally calculated variables")
    setvarlist = []
    for var in varlist:
        dims = (len(bbbv[var])-2)
        var = f"{var}_tmp(xc,yc)={var}_tmp(xc,yc)+{var}(xc,yc)"
        setvarlist.append(var.replace(")", dims*",:"+")"))
    
        
    lines =  [x.replace(" ","; ") for x in wrap(
            " ".join(setvarlist), width=70, 
            break_long_words=False)]
#    lines =  [x.replace(" ","; ") for x in wrap(
#            " ".join([f"{x}_tmp={x}_tmp+{x}" for x in varlist]), width=70, 
#            break_long_words=False)]
    outlines.append(8*" " + "\n        ".join(lines))
    outlines.append(4*" " + "END DO") 
    outlines.append(4*" " + "!$OMP  END PARALLEL DO")

    outlines.append("\n    ! Update global variables")        
    lines =  [x.replace(" ","; ") for x in wrap(
            " ".join([f"{x}={x}_tmp" for x in varlist]), width=70, 
            break_long_words=False)]
    outlines.append(4*" "+"\n    ".join(lines))

    lines =  [x.replace(" ","; ").replace("!", " ") for x in wrap(
            " ".join([f"call!OmpCopyPointer{x}" for x in varlist]), width=70, 
            break_long_words=False)]
    outlines.append(4*" "+"\n    ".join(lines))
    

    outlines.append(f"\n  END SUBROUTINE OMP{subroutine}")
    
    return outlines
    
        



def get_pandf_vars():
    from os import walk
    subroutines = {}
    bbbv = parse_bbbv()
    for root, dirs, files in walk('bbb'):
        for file in files:
            if (("pandf_" in file) and (file[-2:]==".m")) or (file == "convert.m"):
                for subroutine, varlist in get_sourcefile_vars_set(f"bbb/{file}").items(): 
                    subroutines[subroutine] = varlist
    bbbvars = list(bbbv.keys())
    defined_vars = {}
    for subroutine, varlist in subroutines.items():
        defined_vars[subroutine] = list(set(varlist) & set(bbbvars))
    return defined_vars




def parse_bbbv(fname='bbb/bbb.v'):
    data = {}
    begin = False
    with open(fname, 'r') as f:
        for line in f:
            if begin:
                line = line.split("#")[0].strip()
                if ("****" not in line) and (len(line) > 0):
                    var = line.split()[0]
                    if "(" in var:
                        [varname, vardim] = var.split('(')
                        data[varname.lower()] = vardim.replace(")","").split(",")
                    else:
                        data[var] = "int"
                else:
                    continue
            else:
                if "****" in line:
                    begin = True
                else:
                    continue
    return data

def get_sourcefile_vars_set(fname):
    subroutines = {}
    subrutine = None
    locdef = False
    with open(fname) as f:
        for line in f:
            if locdef is True:
                if len(line.strip())>0:
                    if line.strip()[0] == '.':
                        for var in line.replace('.','. ').replace(',', ' ').replace('(', ' ').split():
                            var.strip()
                            if len(var)>0:
                                block['defined'].append(var)
                        continue
            if (line[0].lower() == 'c') or (line[0] == '!') or (line[0] == '*'):
                continue
            line = line.split('#')[0]
            line = line.split('remark')[0]
            line = line.split('xerrab')[0]
            if ('subroutine' in line.lower()) and (line.split()[0].lower() == 'subroutine'):
                subroutine = line.split()[1].split('(')[0].strip()
                subroutines[subroutine] = []
                block = subroutines[subroutine]
            elif 'end subroutine' in line.lower():
                subroutine = None
                block = None
            elif len(line.strip())==0:
                continue
            elif (line.strip()[:4].lower() == 'call'):
                locdef = False
                continue
            elif "=" in line:
                try:
                    if int(line.split("=")[1].strip().replace(".","")) == 0:
                        continue
                except:
                    pass
                line = line.split("=")[0].lower().replace(".","").strip()
                if line[:2] in ['if', 'do']:
                    continue
                elif "(" in line:
                    line = line.split("(")[0]
                subroutines[subroutine].append(line)
            locdef = False
    for subroutine, varlist in subroutines.items():
        subroutines[subroutine] = list(dict.fromkeys(varlist))
    return subroutines




def get_sourcefile_vars_used(fname):
    from uetools import Case
    c=Case()
    subroutines = {}
    subrutine = None
    locdef = False
    with open(fname) as f:
        for line in f:
            if locdef is True:
                if len(line.strip())>0:
                    if line.strip()[0] == '.':
                        for var in line.replace('.','. ').replace(',', ' ').replace('(', ' ').split():
                            var.strip()
                            if len(var)>0:
                                block['defined'].append(var)
                        continue
            if (line[0].lower() == 'c') or (line[0] == '!') or (line[0] == '*'):
                continue
            line = line.split('#')[0]
            line = line.split('remark')[0]
            line = line.split('xerrab')[0]
            if ('subroutine' in line.lower()) and (line.split()[0].lower() == 'subroutine'):
                subroutine = line.split()[1].split('(')[0].strip()
                subroutines[subroutine] = {'packages': [], 'local': [], 'undef': [], 'imported': [], 'defined': [],
                    'missing_packages': [], 'missing_localvars': []} 
                block = subroutines[subroutine]
            elif 'end subroutine' in line.lower():
                subroutine = None
                block = None
            elif len(line.strip())==0:
                continue
            elif (line.strip()[:4].lower() == 'call'):
                locdef = False
                continue
            # Skip package imports
            # Skip declarations and imports
            elif line.split()[0].lower() == 'implicit':
                locdef = False
                continue
            elif line.split()[0].lower() in ['integer', 'real', 'logical', 'external', 'parameter']:
                locdef = True
                for vartype in ['integer', 'real', 'logical', 'external','parameter']:
                    line = line.replace(vartype, ' ')
                for var in line.split(','):
                    var = var.split('(')[0].strip()
                    if len(var)>0:
                        block['defined'].append(var.lower())
                continue
            elif line.split()[0][:9] == 'character':
                line = line.split('*')[1].split()[1:]
                for var in line:
                    block['defined'].append(var.lower())
                
            # Skip package imports
            elif line.strip()[:4].lower() == 'use(':
                block['imported'].append(line.split('(')[1].split(')')[0].lower())
            # Parse lines
            elif subroutine is not None:
                for split in ['.and.', '.or.', '.lt.', '.le.', '.ge.', '.gt.',
                    '=','+','-','*','/', '(',')',',', 'elseif', 'endif','if ','do',
                    '<','>','then',' min',' max', 'else','.eq.', ' log', 'exp',
                    ' mod', ' sqrt', ' abs', '.', ':', 'cos', 'sign', 'enddo',
                    'return', ' end', 'cycle','false', 'true', 'continue',
                    'call', 'goto', 'break']:
                    for substr in [split, ' '+split,
                        split+' ', split.replace('(', ' (')]:
                        line = line.lower().replace(substr,' ')
                for var in line.split():
                    try:
                        float(var)
                        continue
                    except:
                        pass
                    try:
                        block['packages'].append(c.search.get_group(var))
                    except:
                        if (var[0] !='"') and (var[0] !="'"):
                            block['local'].append(var)
            locdef = False

    for blockname, entries in subroutines.items():
        block = subroutines[blockname]
        for varlist in list(block.keys()):
            block[varlist] = list(dict.fromkeys(block[varlist]))
        for var in block['packages']:
            if var.lower() not in block['imported']:
                block['missing_packages'].append(var)
        for var in block['local']:
            if var.lower() not in block['defined']:
                block['missing_localvars'].append(var)
        if len(block['missing_packages']) + len(block['missing_localvars']) > 0:
            block['empty'] = False
        else:
            block['empty'] = True
         
    for subroutine, entries in subroutines.items():
        if entries['empty'] is True:
            continue
        print(f'======= {subroutine} ========')
        if len(entries['missing_packages']) > 0:
            print("Missing imports:")
            for var in entries['missing_packages']:
                print(f"      Use({var})")
        if len(entries['missing_localvars']) > 0:
            print("Missing local vars:")
            print("   " + ", ".join(entries['missing_localvars']) )
            
    for sr, item in subroutines.items():
        print(sr)
 
