import pytest

from berny import Berny, geomlib, optimize
from berny.solvers import MopacSolver

analine_str = """\
14
Aniline
H      1.5205     -0.1372      2.5286
C      0.9575     -0.0905      1.5914
C     -0.4298     -0.1902      1.6060
H     -0.9578     -0.3156      2.5570
C     -1.1520     -0.1316      0.4215
H     -2.2452     -0.2104      0.4492
C     -0.4779      0.0324     -0.7969
N     -1.2191      0.2008     -2.0081
H     -2.0974     -0.2669     -1.9681
H     -0.6944     -0.0913     -2.8025
C      0.9208      0.1292     -0.8109
H      1.4628      0.2560     -1.7555
C      1.6275      0.0685      0.3828
H      2.7196      0.1470      0.3709
"""

cyanogen_str = """\
4

N      3.545830    3.669192    7.228181
C      3.601888    3.624940    6.062501
C      3.671915    3.575700    4.697549
N      3.727670    3.537778    3.532496
"""

ethanol_str = """\
9

C	1.1879	-0.3829	0.0000
C	0.0000	0.5526	0.0000
O	-1.1867	-0.2472	0.0000
H	-1.9237	0.3850	0.0000
H	2.0985	0.2306	0.0000
H	1.1184	-1.0093	0.8869
H	1.1184	-1.0093	-0.8869
H	-0.0227	1.1812	0.8852
H	-0.0227	1.1812	-0.8852
"""

water_str = """\
3

O 0 0 0
H 0 0.7935 0
H 0 0 0.7935
"""


@pytest.fixture
def mopac():
    return MopacSolver()


def ethanol():
    return geomlib.loads(ethanol_str, fmt='xyz'), 5

# TODO: Need to see why this takes 12 steps instead of 11
def aniline():
    return geomlib.loads(analine_str, fmt='xyz'), 12


def cyanogen():
    return geomlib.loads(cyanogen_str, fmt='xyz'), 4


def water():
    return geomlib.loads(water_str, fmt='xyz'), 7


@pytest.mark.parametrize('test_case', [ethanol, aniline, cyanogen, water])
def test_optimize(mopac, test_case):
    geom, n_ref = test_case()
    berny = Berny(geom)
    optimize(berny, mopac)
    assert berny.converged
    assert berny._n == n_ref
