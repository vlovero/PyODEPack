from typing import overload, Callable
from numpy.typing import ArrayLike, Tuple
from numpy import ndarray, inf
from . import methods


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
