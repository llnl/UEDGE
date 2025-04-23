from uetools import Case
c=Case()
subroutines = {}
subrutine = None
locdef = False
with open('bbb/pandf.m') as f:
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
        
