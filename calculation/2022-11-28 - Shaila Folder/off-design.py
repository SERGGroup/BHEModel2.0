# %% ------------ Design parameters ------------------------------------------------------------------->
Y_id=[0.4285451346283881, 0.34819140402298704, 0.2799601412216672, 0.22278180031350817]
m_dot=36.38  # required mass flow to produce 1MW at 22 degrees
#when ambient temperature changes the mass flowrate remain constant/ turbine inlet and outlet condition changes
def fi(rho,p,m=m_dot):
    Fi = m/((p*rho)**0.5)
    return Fi
def beta(Fi,Y):
    beta=1/(1-(Fi**2)*Y)
    return beta

# %% ----------------stage 1------------------------------------------------------------------>

T_in=125
P_in=