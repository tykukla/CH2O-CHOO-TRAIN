#-------------------------------------------------------------------------------#
#             Moist Energy Balance Model (MEBM) -- Phys. Constants              # 
#        Gerard Roe (Univ. Washington) and Tyler Kukla (Colostate Univ.)        #
#-------------------------------------------------------------------------------#


# physical constants  - no choice permitted (but feel free to change to more precise values...)
cp <- 1004          # [J kg-1] heat cap. at const. pressure.
eps <- 0.622        # [unitless] constant in specific humidity relationship (=m_v/m_d)
# Note consider adding T dependency of Lv at some point: Lv = Rv*(T+273.15)^2*a*b/(b+T)^2
Lv <- 2.45e6        # [J kg-1] latent heat of vaporization at 20 degC.
e0 <- 611.2         # [Pa] sat. vap. press. in Pa at 0 degC.
a <- 17.67          # [unitless] constant in the esat formula.
b <- 243.5          # [degC] constant in esat formula.
p0 <- 1.013e5       # [Pa] surface reference pressure.
g <- 9.81           # [m s-2] Accel due to gravity
Rv <- 461           # [J kg-1] gas constant for water vapor. Use for CC relationship.
Re <- 6.37e6        # [m] Earth's radius.
rho <- 1e3          # [kg m-3] density of water
