import numpy as np
import EA2D

def ising_run(nb_steps, T, H, J):
    """
    Utility python wrapper around Fortran function EA2D for fast computations
    involving the Ising model.

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
    M : np.array, float
      energy and magnetisation at each step

    Pedro Mediano, Mar 2020
    """

    # Dummy array for return values
    V = np.zeros(2*nb_steps)

    # Run simulation
    EA2D.ea2d(V, T/10, H, J) # modified it to T/10, similar to the original script

    # Get returned values and exit
    E = V[:nb_steps]
    M = V[nb_steps:]

    return E, M


if __name__ == '__main__':
    nb_steps = 100000
    for size in [11, 21, 51, 101, 121, 151]:
        for T in [2, 6, 8, 9, 10, 11, 12, 15, 20, 30, 40]:
            print("Iteracion size:{}, T:{}".format(size,T))
            H = np.asfortranarray(np.loadtxt('./IsingParams/N{}H.txt'.format(size)))
            J = np.asfortranarray(np.loadtxt('./IsingParams/N{}J.txt'.format(size)))

            E, M = ising_run(nb_steps, T, H, J)
            np.savetxt("Results/T{}N{}E.txt".format(T,size),E) #Save the energy 
            np.savetxt("Results/T{}N{}M.txt".format(T,size),M) #Save the Magnetization
    print('Done!')

