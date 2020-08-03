#====================================================
# tools
def summary(updates = [], prefix = 'datas'):
    datas = {}
    calc = Espresso()
    print('Reading.....')
    for update in updates:
        cwd = os.getcwd()
        for i,j,y in os.walk(update):
            output = is_espresso(i)
            if output:
                os.chdir(i)
                print('dire:', i)
                calc.directory = cwd + '/' + i
                calc.prefix = output[0:-4]
                try:
                    calc.results = {}
                    calc.read_results()
                    t = calc.get_time()
                    calc.results['time'] = t
                    datas[i] = calc.results
                except Exception as e:
                    print('='*30, '\n', i, e)
            os.chdir(cwd)
    with open('%s.pickle' % prefix, 'wb') as f:
        pickle.dump(datas, f)
    print('Finished')

def is_espresso(path):
    '''
    check espresso 
    '''
    dirs = os.listdir(path)
    # print(dirs)
    # flag = True
    for qefile in ['.pwo']:
        flag = False
        for file in dirs:
            if qefile in file:
                return file
        if not flag:
            return False
    # return flag
def grep_valence_configuration(pseudopotential):
    """
    Given a UPF pseudopotential file, find the valence configuration.
    
    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    3S  1  0  2.00      0.700      1.200    -6.910117

    """
    orbitals = ['S', 'P', 'D', 'F']
    valence = {}
    with open(pseudopotential) as psfile:
        lines = psfile.readlines()
        for i in range(len(lines)):
            if 'valence configuration:' in lines[i].lower():
                j = i + 2
                ob = lines[j].split()[0]
                while ob[1] in orbitals:
                    valence[ob] = lines[j].split()[3]
                    j += 1
                    ob = lines[j].split()[0]
                return valence
    if not valence:
        raise ValueError('Valence configuration missing in {}'.format(pseudopotential))

