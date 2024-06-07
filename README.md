# Quercetin/DMF Surface Analysis

Supporting code for <placeholder for DOI>.

---

# Using Script

The `quercetin_DMF_surface_analysis.py` script can be used through the command line.

First we activate the miniconda that contains the `csd-python-api`.

```commandline
"C:\Program Files\CCDC\ccdc-software\csd-python-api\miniconda\Scripts\activate.bat"
```

The script can then be run for on a specific surface of a given refcode from the CSD or a `.cif` file by specifying the
surface of interest as three integers h, k, and l:

```commandline
python quercetin_DMF_surface_analysis.py REF_CODE|file_name.cif h k l -o 0.00
```

Output will then be printed to the console.

By default, the script will report values for an offset of 0.00 Å from the specified surface. This offset can be changed
using the `-o` command followed by a float. For example, the following command:

```commandline
python quercetin_DMF_surface_analysis.py HXACAN 2 0 0 -o -2.951 
```

will print the following output:

```commandline
O.2 : 6.36
O.3 : 18.958
H : 54.793
N.am : 1.034
C.3 : 12.862
C.ar : 42.315
{'O.2': 6.36, 'O.3': 18.958, 'H': 54.793, 'N.am': 1.034, 'C.3': 12.862, 'C.ar': 42.315}
```

# Dependencies

- [CSD-Python-API](https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd-python-api/)

# License Requirements

- CSD-Particle

_CSD-Particle is available by default as part of an academic CSD-Enterprise license._