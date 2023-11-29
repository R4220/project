
def setup():
    '''
    This function extracts the setup information from Setup.txt

    Args:
        None

    Returns:
        filename: the name of the pwo file that we'll convert into a xyz file
    '''

    with open('Setup.txt', 'r') as fset:
        for line in fset:
            if 'Filename:' in line:
                filename = line.split()[1]

    return filename

def xyz_lines(fin):
    for idx,line in enumerate(fin):
        fout.writelines()
    return

if __name__ == "__main__":
    filename = setup()
    with open(filename + '.xyz', "w+") as fout:
        with open(filename + '.pwo', 'r') as fin:
            xyz_lines(fin)
        
        #usa le variabili ad errore, magari con un dizionario, a cui soino collegate
        #fout.writelines(["%s\n" % i for i in gen_list_for_xyz(file_to_open + '.pwo')])

 #   '''file_to_open = str(input('Insert filename: ____.pwo \n'))'''