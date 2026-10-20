"""Independent P1 generalized eigenproblem and alternate coupled root solve."""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import root
from scipy.special import ai_zeros
from quantum_mos import HBAR,M0,Q


def finite_element_spectrum(z,potential_eV,mass,states=6):
    """Consistent mass matrix; exact Gauss integration of linear V times Ni Nj."""
    z=np.asarray(z); potential_eV=np.asarray(potential_eV)
    n=len(z); k=np.zeros((n,n)); m=np.zeros((n,n))
    kinetic=HBAR**2/(2*mass*M0*Q)
    gauss=np.array([.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)])
    for i in range(n-1):
        h=z[i+1]-z[i]; ids=np.array([i,i+1])
        local_k=kinetic/h*np.array([[1.,-1.],[-1.,1.]])
        local_m=h/6*np.array([[2.,1.],[1.,2.]])
        for t in gauss:
            shape=np.array([1-t,t]); v=shape@potential_eV[ids]
            local_k+=h/2*v*np.outer(shape,shape)
        k[np.ix_(ids,ids)]+=local_k; m[np.ix_(ids,ids)]+=local_m
    e,w=eigh(k[1:-1,1:-1],m[1:-1,1:-1],subset_by_index=(0,states-1))
    full=np.zeros((n,states)); full[1:-1]=w
    return e,full


def box_levels(length_m,mass,states=4):
    j=np.arange(1,states+1)
    return HBAR**2*np.pi**2*j*j/(2*mass*M0*Q*length_m**2)


def triangular_levels(slope_eV_nm,mass,states=4):
    c=HBAR**2/(2*mass*M0*Q)*1e18
    return -ai_zeros(states)[0]*(c*slope_eV_nm**2)**(1/3)


def coupled_root(model,vg):
    """Same discrete equations, different nonlinear algorithm; not a new physics model."""
    classical=model.solve(vg,"classical")
    evaluations=0
    def residual(psi):
        nonlocal evaluations
        evaluations+=1
        n,_,_=model.quantum_density(psi)
        return model.band_solve(model.residual(psi,n,vg))
    solved=root(residual,classical["potential_V"],method="hybr",options={"xtol":1e-10})
    error=float(np.max(abs(residual(solved.x))))
    if not solved.success or error>1e-8:
        raise RuntimeError(f"Alternate coupled root failed: {solved.message}; correction={error}")
    n,_,_=model.quantum_density(solved.x)
    return dict(potential_V=solved.x,Ns_m2=float(model.weights@n),correction_V=error,evaluations=evaluations)
