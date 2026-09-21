<!-- Title -->
<h1 align="center">
  Quantizing Magnetic Fields in Magnetars (QMFM)
</h1>

<!-- Description -->
QMFM is a finite volume code for studying magnetothermal evolution in neutron star crusts, and includes the effects of Landau quantization on thermodynamic quantities and transport coefficients. It is written in C++ and parallelized with MPI.

<table style="background-color:#FFFFFF;">
  <tr>
    <th width="50%">
        <figure>
          <img src="https://raw.githubusercontent.com/QMFM/extra/BannerMagTherm.png">
          <figcaption>Magnetothermal evolution</figcaption>
        </figure>
      </a>
    </th>
    <th width="50%">
        <figure>
          <img src="https://raw.githubusercontent.com/QMFM/extra/BannerCurrent.png">
          <figcaption>Current density</figcaption>
        </figure>
      </a>
    </th>
  </tr>
  <tr>
</table>

#Installation instructions and requirements

QMFM requires the following external libraries:

* GNU Scientific Library: <https://www.gnu.org/software/gsl/>
* FFTW3: <https://www.fftw.org/>
* SHTns: <https://nschaeff.bitbucket.io/shtns/>
* Boost: <https://www.boost.org/>
* OpenMPI: <https://www.open-mpi.org/software/ompi/v5.0/>
* HDF5: <https://github.com/HDFGroup/hdf5/releases>

To install QMFM, in terminal run

```terminal
git clone https://github.com/PRauNS/QMFM.git
mkdir build-dir
cd build-dir
cmake ..
make
```

which will clone the repo and generate the executable "main" in the folder "QMFM". An example CMake file (CMakeLists.txt) is included; it may need to be adjusted based on where the above libraries are installed on the user's system.

#Running the code

After installation, the code is run using
```terminal
./main
```
from the QMFM directory. The options for running simulations are specified within the SimSetup.in file; the options are described by comments in this file. Note that adding additional entries to this file that are not preceded by hashmarks, indicating a comment, will likely cause the code to crash.

The only compatible equation of state options at present are the BSk24 (<https://doi.org/10.1093/mnras/stz800>) and SLy4 (<https://www.aanda.org/articles/aa/abs/2001/46/aa1755/aa1755.html>); data files for a 1.4 solar mass model for each are included in the repo. These are loaded using the options

CrustEOS: TOVBSk24Crust.dat
CoreEOS: TOVBSk24Core.dat

or 

CrustEOS: TOVSLy4Crust.dat
CoreEOS: TOVSLy4Core.dat

respectively.

The initial magnetic field at present only of the form described in; the strength of its poloidal and toroidal components is specified (in G) by options (e.g., for a 10^14 G poloidal field at the poles and a 10^14 G maximum toroidal field)

B_pol_init: 1e14
B_tor_init: 1e14

#Developers

Peter B. Rau <https://github.com/PRauNS>
