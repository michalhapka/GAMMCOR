# GammCor plugin for Quantum Package
Authors: Michal Hapka, Anthony Scemama, Kasia Pernal

This is a plugin for the Quantum Package program: 
https://quantum-package.readthedocs.io/en/master/index.html

The plugin interfaces Quantum package and the GammCor code:
https://qchem.gitlab.io/gammcor-manual/index.html

via the TREXIO library:
https://github.com/trex-coe/trexio

## Installation
#### 1. Clone the repository to your Quantum Package /plugins directory
```
git clone git@gitlab.com:michalhapka/gammcor_qp2_plugin.git
```

#### 2. Source the quantum_package.rc file

#### 3. Install the plugin
```
qp_plugins install gammcor_plugin
```

## Usage
```
qp set_file <your_EZFIO_dir>
qp set gammcor_plugin cholesky_tolerance 1.e-5
qp set gammcor_plugin trexio_file <HDF5_filename>
qp run export_gammcor >  output.out
qp run gammcor_plugin >> output.out
```

## License
GNU GENERAL PUBLIC LICENSE, Version 2, June 1991
