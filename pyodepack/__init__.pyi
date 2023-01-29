from typing import overload, Callable
from numpy.typing import ArrayLike, Tuple
from numpy import ndarray, inf
from . import methods


class PyCapsule:
    pass


@overload
def solveIVP(f: Callable, y0: ArrayLike, teval: ArrayLike, args: tuple = None, out: ndarray = None,
             isvoid: bool = False, atol: float = 1.49012e-8, rtol: float = 1.49012e-8, min_step: float = 1e-6,
             max_step: float = inf, window: ArrayLike = None, method: int = methods.LSODA, full_output: bool = False) -> Tuple[ndarray, int]:
    """Integrate a system of ODEs.

    Parameters
    ----------
    f : Callable
        RHS of system
    y0 : ArrayLike
        _Initial conditions
    teval : ArrayLike
        points to evaluate solution at (must include initial time point t0)
    args : tuple, optional
        extra arguments for f(t, y), by default None
    out : ndarray, optional
        preallocated output array, by default None
    isvoid : bool, optional
        f does not return anything (performance booster), by default False
    atol : float, optional
        absolute error tolerance, by default 1.49012e-8
    rtol : float, optional
        relative error tolerance, by default 1.49012e-8
    min_step : float, optional
        minimum allowable step, by default 1e-6
    max_step : float, optional
        maximum allowable step, by default inf
    window : ArrayLike, optional
        computational window solution must stay inside, by default None
    method : int, optional
        which integration method to use, by default methods.LSODA ( == 0)
    full_output : bool, optional
        return dict with integration stats, by default False

    Returns
    -------
    Tuple[ndarray, int]
        ndarray with shape (len(teval), len(y0)) containing solution
        and number of successful check points (from teval) reached.
    Tuple[ndarray, int, Dict[str, Any]]
        If full_output == True, also return dictionary with integtation stats.
    """
    ...


@overload
def ODEStep(f: Callable, y0: ArrayLike, t0: float, tf: float, args: tuple = None, out: ndarray = None,
            isvoid: bool = False, atol: float = 1.49012e-8, rtol: float = 1.49012e-8, min_step: float = 0.0,
            max_step: float = inf, window: ArrayLike = None, method: int = methods.RKV65E, full_output: bool = False,
            h0: float = 0.0, first_step: bool = True, work: ndarray = None, stepper: PyCapsule = None) -> ndarray:
    """"""

    """Integrate a system of ODEs from t0 to tf.

    Parameters
    ----------
    f : Callable
        RHS of system
    y0 : ArrayLike
        _Initial conditions
    t0 : float
        starting time
    tf : float
        finishing time
    args : tuple, optional
        extra arguments for f(t, y), by default None
    out : ndarray, optional
        preallocated output array, by default None
    isvoid : bool, optional
        f does not return anything (performance booster), by default False
    atol : float, optional
        absolute error tolerance, by default 1.49012e-8
    rtol : float, optional
        relative error tolerance, by default 1.49012e-8
    min_step : float, optional
        minimum allowable step, by default 1e-6
    max_step : float, optional
        maximum allowable step, by default inf
    window : ArrayLike, optional
        computational window solution must stay inside, by default None
    method : int, optional
        which integration method to use, by default methods.LSODA ( == 0)
    full_output : bool, optional
        return dict with integration stats, by default False
    h0 : float, optional
        initial time step, by default 0.0
    first_step : bool, optional
        first step being taken by stepper, by default True
    work : ndarray, optional
        work array for method, by default None
    stepper : PyCapsule, optional
        PyCapsule containing pointer to C++ stepper, returned after first step, by default None

    Returns
    -------
    ndarray
        ndarray with shape (len(y0),) containing solution at tf.
    Tuple[ndarray, Dict[str, Any]]
        If full_output == True, also return dictionary with integtation stats.
    """
    ...


@overload
def parse_stepper(stepper: PyCapsule):
    """Convert C++ stepper object into dictionary of useful information.

    Parameters
    ----------
    stepper : PyCapsule
        PyCapsule containing pointer to stepper.
    """
    ...


@overload
def parse_work(work: ndarray, method: int, node: int) -> dict:
    """Convert work array from full_output solve into readable dict.

    Parameters
    ----------
    work : ndarray
        work array from ODEStep
    method : int
        method from `methods` module
    node : int
        number of odes of system being solved

    Returns
    -------
    dict
        dictionary of parsed work
    """
    pass
