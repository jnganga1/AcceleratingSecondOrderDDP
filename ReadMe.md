Requirements
------------

-   This code has a dependency on CASADI. The documentation of casadi can be
    found here:* https://web.casadi.org* We include both windows casadi
    implementation in *"casadi-windows-matlabR2016a-v3.5.3"* folder and Linux
    implementation  in *“casadi-linux-matlabR2014b-v3.5.5”* folder. To use them
    within matlab, they need to be added to PATH.

-   This code also uses Spatial_V2, *http://royfeatherstone.org/spatial/v2/*, to
    create models -\> this is included in "Spatial_v2" folder.

-   To interface spatial_v2 with casadi, we implement some special  Spatial_v2
    functions using casadi functionality as in *"spatial_v2_casadi"* folder

Folder Organization
-------------------

-   *[Algorithm]*** **- contains DDP algorithms and all methods of attaining
    derivatives

-   *[casadi-windows-matlabR2016a-v3.5.3]* - contains casadi's Windows
    implementation

-   *[casadi-linux-matlabR2014b-v3.5.5]* - contains casadi's Linux
    implementation

-   *[DDP_Regular_ABA_Methods - LBRModel]*** **- implements DDP for a 7-link
    KUKA LBR

-   *[Spatial_v2]* - Cointains spatial v2 implementation.

-   *[Spatial_V2_casadi]* - Implements certain spatial v2 function using the
    casadi toolkit

Important Files
---------------

-   *[Tester_partialsEval.m]* - evaluates a single evaluation of dynamics
    partials for varying number of DoFs. Produces Fig. 4 (Fig.1 is a s subset of
    Fig. 4)

-   *[Tester_VaryLinks.m]* - evaluates variation of links vs time. Produces Fig.
    1 and Fig. 5

-   *[Tester_VaryingInitialConditions.m]* - randomizes the initial condition for
    n = 7. Produces Fig. 6

-   *[Tester_randomizeLBR.m]* - randomizes controller for LBR KUKA robot and
    plots its results to produce Fig. 7 and Fig 8
