def get_subroutine_vars(subroutinelist):
    defined_vars = get_pandf_vars()
    var = []
    for subroutine in subroutinelist:
        var = var + defined_vars[subroutine]
    return var

def get_subpandf1():
    from uetools import Case
    getpkg = Case().search.get_group
    pkgs = {}
    varlist = get_subroutine_vars(['convsr_vo', 'convsr_aux'])
    for var in varlist:
        pkg = getpkg(var)
        if pkg not in pkgs:
            pkgs[pkg] = []
        pkgs[pkg].append(var)

    bbbv = parse_bbbv()
    vardef = []
    for var in varlist:
        vardef.append("{}_tmp({})".format(var, ",".join(bbbv[var])))

    for pkg, var in pkgs.items():
        print("    USE {}, ONLY: {}".format(pkg, ", &\n    &    ".join(var)))


    print("\n\n\n! Define local variables")
    print("real:: {}".format(", &\n    &      ".join(vardef)))

    print("\n\n\n! Initialize arrays to zero")        
    for var in varlist:
        print(f"{var}_tmp = 0.")

    print("\n\n\n          !$omp &   REDUCTION(+:yldottot,{})".format(\
            ', &\n          !$omp       & '.join([f"{x}_tmp" for x in varlist])))        

    print("\n\n\n! Update locally calculated variables")            
    for var in varlist:
        print(f"{var}_tmp = {var}")

    print("\n\n\n! Update global variables")        
    for var in varlist:
        print(f"{var} = {var}_tmp")


    
        



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
