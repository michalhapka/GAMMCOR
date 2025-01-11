# PYSCF - GAMMCOR Interface Description

PYSCF works in INCORE or simple THC mode. For THC, you need to set Cholesky OTF in input.

## 1. Setup and Example

**1.1Example Input (eth.py)**

```python
#!/usr/bin/env python
import sys
sys.path.append('/path/to/your/gammcor/')  # Replace with your path
import numpy as np
from pyscf import gto, scf, mcscf
from gammcor.interface.pyscf_functions import *

# Create molecule
mol = gto.M()
mol.atom = '''
C   0.0000   0.0000   0.6695 
C   0.0000   0.0000  -0.6695 
H   0.0000   0.9289   1.2321 
H   0.0000  -0.9289   1.2321 
H   0.0000   0.9289  -1.2321 
H   0.0000  -0.9289  -1.2321
'''
mol.basis = {'C': 'def2-tzvp', 'H': 'def2-tzvp'}
mol.symmetry = True
mol.build()

# Run HF
myhf = scf.RHF(mol)
myhf.kernel()

# Run CASSCF
norb = 2  # active orbitals
nelec = 2 # active electrons
mycas = mcscf.CASSCF(myhf, norb, nelec)
mycas.natorb = True  # use natural orbitals
mycas.kernel()

# Dump files for GAMMCOR
get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)
```

**1.2 Required Setup**

Before running any calculations:

1. Replace this with your path to python script:
```python
sys.path.append('/path/to/your/gammcor/')
```

2. Replace this with  the correct path to basis set files:
```python
basis_file = '/path/to/your/gammcor/basis-for-pyscf/basis.nw'
```

### 2 Running Calculations

**2.1 Basic CASSCF Setup:**

- Always use `natorb=True` for single-state calculations
- Canonicalization is on by default (recommended)
- For State-Averaged calculations naturalization is done in Gammcor

**2.2 Dumping Files for GAMMCOR:**

```python
# After CASSCF calculation, add:
get_data_for_gammcor(mol, myhf, mycas, dump_eri=True)
```
   - Set `dump_eri=True` if you want to dump two-electron integrals
   - Use `dump_eri=False` or simply  `get_data_for_gammcor(mol, myhf, mycas)` for direct (non-incore) calculations

## 3 File Reading Logic

**3.1 File Naming Patterns**

1. Single-state calculations use simple names (e.g., `auxdata.txt`, `rdm2_aaaa.bin`)
2. State-averaged calculations use numbered suffixes in Molpro style (e.g., `auxdata_1.1.txt`, `rdm2_aaaa_1.1.bin`)




**3.2 Binary Files**

- `CAONO.bin`: Orbital coefficients matrix (float64)
  
  - When `natural = 0`: Contains MO coefficients
  - When `natural = 1`: Contains NO coefficients
  
- `HCore.bin`: One-electron integrals matrix (float64)
  
  - In AO basis
  
- `rdm2_aaaa.bin`: Alpha-Alpha reduced density matrix (active space)

- `rdm2_abab.bin`: Alpha-Beta reduced density matrix (active space)

- `rdm2_full.bin`: Full 2-RDM (2.0 * (dm_aaaa + dm_abab))
  
- all rdm.bin files
  
  - When `natural = 0`: In MO basis
  - When `natural = 1`: In NO basis
  
- `TWOEl.bin`: Two-electron integrals (optional, if dump_eri=True)
  
  - In AO basis
  
    

**3.3 Text Files**

- `auxdata.txt`: Contains basic calculation parameters (one per line):
  1. Number of basis functions (nbasis)
  2. Number of inactive orbitals (NI)
  3. Number of active orbitals (NA)
  4. Number of virtual orbitals (NV)
  5. CASSCF total energy
  6. Nuclear repulsion energy
  7. Number of electrons
  8. Natural orbitals flag (1 if used, 0 if not)

- `rdm2.dat`: Full 2-RDM in text format (always in NO basis)

**3.4 Molden Files**

- `moldenhf.inp`: Hartree-Fock orbitals in AO basis
- `molden_mcscf.inp`: CASSCF orbitals

### Reading Cases - this is done automatically. Below is simply the description.

1. **No State Specified in Gammcor input** (when `nst(1) < 0`):
   - If only `auxdata.txt` exists: Uses it, assumes state 1.1
   - If only `auxdata_x.y.txt` exists: Uses it and set nst to x.y
   - If `auxdata.txt` and `auxdata_x.y.txt` files exist: Uses `auxdata.txt`, assumes state 1.1
   - If multiple `auxdata_x.y.txt` files exist and no `auxdata.txt`: Exits with error.

2. **State Specified in Gammcor input** (when `nst(1) > 0`):
   - Looks for either `auxdata.txt` or `auxdata_x.y.txt` matching the specified state
   - Exits if both or neither exist
   - Uses the specified state numbering for reading other files

3. **State Number Setup (anst)**:

When no state is specified (nst(1) < 0):

   - Single `auxdata.txt`: Sets anst = [1,1]
   - Single state-specific file (e.g., `auxdata_2.3.txt`): Sets anst to those numbers [2,3]
   - Multiple state files: Errors unless `auxdata.txt` exists (then uses [1,1])

When state is specified in input to x.y (nst(1) > 0):

   - Sets anst to the specified state numbers
   - Must find exactly one matching file (either `auxdata.txt` or `auxdata_x.y.txt`)
