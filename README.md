# PyODEPack
Fast python wrappers for ODE solvers written in C++

## Available Solvers
 * LSODA
 * Verner 6(5) (efficient)
 * Bogacki-Shampine (RK23)
 * Verner 9(8) (efficient)
 * RadauIIA-5
 * LOVERO competitor to LSODE (nonstiff)

## Requirements
 * setuptools
 * numpy
 * scipy (needed to load LAPACK routines at runtime)

## Install
`SYSTEM_VERSION_COMPAT=1 python -m pip install git+https://github.com/vlovero/PyODEPack.git`
The `SYSTEM_VERSION_COMPAT` environment variable is to fix the
MacOS issue that causes clang to not like the `march=native` compiler flag.
 
 ## Example Usage
 ```Python
 >>> import numpy as np
>>> from pyodepack import solveIVP, methods
>>> from pprint import pprint
>>>
>>> def f(t, z):
...     return np.array([z[1], -z[0]])
...
>>> 
>>> def g(t, z, out, args):
...     k1, k2 = args
...     out[0] = -k1**2 * z[0]
...     out[1] = -k2**2 * z[1]
...
>>> y, _, info = solveIVP(f, [0, 1], [0, 20], method=methods.RKV65E, full_output=True)
>>>
>>> pprint(info)
{'fails': 0,
 'feval': 683,
 'hout': 0.4224769620224126,
 'jeval': 0,
 'message': 'success',
 'numlu': 0,
 'steps': 85}
>>> np.abs(y[:, 0] - np.sin([0, 20]))
array([0.0000000e+00, 3.5364417e-09])
>>> y, _, info = solveIVP(g, [0, 1], [0, 20], args=(1, 20), isvoid=True, full_output=True)
>>>
>>> pprint(info)
{'fails': 0,
 'feval': 213,
 'hout': 79.76670760332844,
 'jeval': 5,
 'message': 'success',
 'numlu': 0,
 'steps': 104}
 ```