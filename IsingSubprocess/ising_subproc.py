"""
Subprocess interface to Ising Metropolis Fortran code.

This interface is designed to keep as much of the original Fortran code intact.
To this end, it follows the same proceduce of changing variable values in the
source and recompiling. Then, it calls the executable via a subprocess and
captures stdout.


Pedro Mediano, Apr 2020
"""

import numpy as np
import subprocess as sp
import tempfile
from io import StringIO
from os.path import isfile


def ising_subproc(nb_steps, T, H, J):
    """
    Utility python wrapper around Fortran function EA2D for fast computations
    involving the Ising model.

    WARNING: There cannot be any files called 'variables.mod' or 'mt95f4.mod'
    in the current directory. If there are, gfortran will use the wrong ones
    and the code will not work as expected.

    Parameters
    ----------
    nb_steps : int
      number of steps to run the simulation for
    T : float
      temperature
    H : np.array, float
      external field array
    J : np.array, float
      coupling matrix

    Returns
    -------
    E : np.array, float
      energy at each step

    Pedro Mediano, Apr 2020
    """

    ## Parameter checks
    n1, n2 = J.shape
    n3 = len(H)
    if n1 != n2 or n1 != n3 or n2 != n3:
        raise ValueError('Arrays H and J must have consistent sizes.')
    if isfile('variables.mod') or isfile('mt95f4.mod'):
        raise RuntimeError("There cannot be any files called 'variables.mod'"+
                " or 'mt95f4.mod' in the current folder. Please rename or remove them.")

    ## Write H and J matrices to a file
    hFile = tempfile.NamedTemporaryFile(delete=True)
    np.savetxt(hFile, H, fmt='%f')
    hFile.file.flush()

    jFile = tempfile.NamedTemporaryFile(delete=True)
    np.savetxt(jFile, J, fmt='%f')
    jFile.file.flush()
    
    ## Load original source code
    with open('specific_heat.f90', 'r') as f:
        src = f.read()

    ## Search-replace the code to pass arguments (ugh)
    src = src.replace('100000', str(nb_steps))
    src = src.replace('Largo=101', 'Largo=%i'%len(H))
    src = src.replace('T=8', 'T=%i'%(10*T))
    src = src.replace('HFILENAME', hFile.name)
    src = src.replace('JFILENAME', jFile.name)

    d = tempfile.TemporaryDirectory(prefix='./')

    ## Write modified source
    with open(d.name+'/RX.f90', 'w') as f:
        f.write(src)

    ## Compile code with new configuration
    p = sp.Popen(['gfortran', '-O3', '-o'+d.name+'/a.out', '-J'+d.name, d.name+'/RX.f90'])
    p.wait()
    if p.returncode != 0:
        raise RuntimeError('Fortran code failed to compile')

    ## Run compiled code
    p = sp.Popen([d.name+'/a.out'], stdout=sp.PIPE)
    out, err = p.communicate()
    E = np.loadtxt(StringIO(out.decode()))
    
    ## Delete temporary files and return
    hFile.close()
    jFile.close()
    d.cleanup()

    return E


if __name__ == '__main__':
    nb_steps = 100000
    for size in [11, 21, 51, 101, 121, 151]:
        for T in [2, 6, 8, 9, 10, 11, 12, 15, 20, 30, 40]:
            print("Iteracion size:{}, T:{}".format(size,T))
            H = np.loadtxt('./IsingParams/N{}H.txt'.format(size))
            J = np.loadtxt('./IsingParams/N{}J.txt'.format(size))
            E = ising_subproc(nb_steps, T*0.1, H, J)
            np.savetxt("Results/T{}N{}E.txt".format(T,size),E) #Save the energy
    print('Done!')

