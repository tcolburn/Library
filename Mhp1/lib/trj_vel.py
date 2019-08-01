# /opt/schrodinger/run python

>>> from schrodinger.application.desmond import framesettools
>>> fs = framesettools.FrameSet('run.dtr')
>>> f0=fs[0]
>>> print f0.__labels__

    ['BUILDCLASS', 'CHEMICALTIME', 'CREATOR', 'ELAPSED', 'ENERGY', 'EX_ENERGY',
    'FORCE_ENERGY', 'FORMAT', 'GID', 'HOME_BOX', 'ISROGUE', 'KERNEL',
    'KIN_ENERGY','NX', 'NY', 'NZ', 'POSN', 'POT_ENERGY', 'PRESSURE',
    'PRESSURETENSOR', 'PROCESSOR', 'PROVENANCE', 'TEMPERATURE',
    'TEMPERATURE_PER_GROUP', 'TITLE', 'TOPOLOGY', 'VELOCITY', 'VERSION',
    'VOLUME']

>>> for frame in fs:
...     v0 = frame.VELOCITY[0:3]
...     print '   ', frame.CHEMICALTIME,'atom 0 has velocity', v0
...     v1 = frame.VELOCITY[3:6]
...     print '   ', frame.CHEMICALTIME,'atom 1 has velocity', v1
...

    [ 240.] atom 0 has velocity [-5.90974188 -3.00732279  3.18637037]
    [ 240.] atom 1 has velocity [-6.66847563 -0.52196407  7.14664888]
    [ 480.] atom 0 has velocity [ 2.80623078  3.28876615  2.06376624]
    [ 480.] atom 1 has velocity [-19.10564041  -2.46500659  10.8352108 ]
    [ 720.] atom 0 has velocity [-10.43045998   0.31181553   0.17032397]
    [ 720.] atom 1 has velocity [-11.26202393 -11.438097    21.30398178]
    [ 960.] atom 0 has velocity [-0.34902263  0.33109778 -1.17421675]
    [ 960.] atom 1 has velocity [  2.74298549 -20.17019463 -14.64497852]
    [ 1200.] atom 0 has velocity [-5.85332203 -4.33602047  4.13961029]
    [ 1200.] atom 1 has velocity [ 17.52938271   3.54563594   7.55537558]
    [ 1440.] atom 0 has velocity [ 2.38425136 -2.21108794 -6.22821856]
    [ 1440.] atom 1 has velocity [ 36.20913696  26.73501205  -6.73337841]

>>> for frame in fs:
...     v0 = frame.POSITION[0:3]
...     print '   ', frame.CHEMICALTIME,'atom 0 has position', v0
...     v1 = frame.POSITION[3:6]
...     print '   ', frame.CHEMICALTIME,'atom 1 has position', v1
...
    [ 240.] atom 0 has position [  3.57056475  12.93109131   8.99171734]
    [ 240.] atom 1 has position [  4.54565811  12.73975372   9.29860878]
    [ 480.] atom 0 has position [  1.11764598  15.07661438   5.73010445]
    [ 480.] atom 1 has position [  1.34935832  15.55605316   6.62343931]
    [ 720.] atom 0 has position [  3.11503911  12.81669331   4.60396338]
    [ 720.] atom 1 has position [  3.1732831   13.72325134   5.11028385]
    [ 960.] atom 0 has position [  3.55293894  17.15921593   3.49656224]
    [ 960.] atom 1 has position [  4.41087914  16.93230629   4.03882647]
    [ 1200.] atom 0 has position [  6.85661316  18.61204529  -5.74965858]
    [ 1200.] atom 1 has position [  6.53081751  19.59861565  -5.7958622 ]
    [ 1440.] atom 0 has position [ 10.465312    18.97945404  -8.14519501]
    [ 1440.] atom 1 has position [ 10.66126251  18.73317909  -9.13643265]