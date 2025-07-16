![q-e-logo](logo.jpg)

This is the Optimal Basis Function (OBF) code based on Quantum Espresso. The code is used as a post-processing tool to generate wavefunction files in QE format (wfc#.dat) from minimal k-points in either scf/nscf calculation from Quantum Espresso. This code utilizes standard QE module to read/write wavefunction files.

The OBF source code is located in OBF directory.

## USAGE
The user can use the OBF software by compiling
```
./configure [options]
make obf
``` 

Successfully compiling the OBF code will generate 3 binaries: obf_basis.x, obf_ham.x and obf_diag.x which are used to build optimal basis function, arbitary k-dependent Hamiltonian and finally QE wavefunction files wfc.dat.
## PACKAGES

obf_basis.x : Builds optimal basis function from input QE scf/nscf calculation.

obf_ham.x : Builds k-dependent Hamiltonian using optimal basis function.

obf_diag.x : Diagonalizes the k-dependent Hamiltonian and constructs wfc.dat

