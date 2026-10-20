"""Equilibrium electron-only (001) thin-Si double-gate Schrodinger-Poisson.

SI internally except energies in eV and input lengths in nm. Infinite walls;
fixed, fully ionized acceptors; no holes or exchange-correlation. EF=0.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded
from scipy.special import expit
from numpy.polynomial.legendre import leggauss

Q=1.602176634e-19
KB=1.380649e-23
HBAR=1.054571817e-34
M0=9.1093837139e-31
EPS0=8.8541878128e-12
EPS=11.7*EPS0
# Valley degeneracy excludes spin. 2D occupation below already includes spin=2.
VALLEYS=(("Delta2",2,0.916,0.190),("Delta4",4,0.190,np.sqrt(0.916*0.190)))
GX,GW=leggauss(384)


@dataclass(frozen=True)
class Device:
    silicon_nm: float=10.0
    oxide_nm: float=2.0
    temperature_K: float=300.0
    acceptors_cm3: float=1e16
    ec_reference_eV: float=0.35
    gate_offset_V: float=0.0

    def __post_init__(self):
        if not np.all(np.isfinite(list(self.__dict__.values()))):
            raise ValueError("Device parameters must be finite")
        if not 3<=self.silicon_nm<=20 or not .7<=self.oxide_nm<=6:
            raise ValueError("Teaching range: silicon=3..20 nm, oxide=.7..6 nm")
        if not 150<=self.temperature_K<=400 or not 0<=self.acceptors_cm3<=1e17:
            raise ValueError("Teaching range: T=150..400 K, NA=0..1e17 cm^-3")
        if not .15<=self.ec_reference_eV<=.7 or abs(self.gate_offset_V)>1:
            raise ValueError("Use Ec0=.15.. .7 eV and |gate offset|<=1 V")

    @property
    def vt(self): return KB*self.temperature_K/Q

    @property
    def cox(self): return 3.9*EPS0/(self.oxide_nm*1e-9)

    @property
    def na(self): return self.acceptors_cm3*1e6


def fermi_half(eta):
    """F_1/2 and its eta derivative by transformed 384-point quadrature.

    t=x^2 removes the sqrt endpoint. Tail cutoff is 50 kT above max eta.
    The same Fermi-Dirac convention is used for all classical valley groups.
    """
    eta=np.asarray(eta,dtype=float)
    if not np.all(np.isfinite(eta)) or np.max(eta)>80 or np.min(eta)<-100:
        raise ValueError("Fermi integral outside supported eta range [-100,80]")
    upper=np.sqrt(max(float(np.max(eta)),0)+50)
    x=(GX+1)*upper/2; weights=GW*upper/2
    occupation=expit(eta[...,None]-x*x)
    scale=4/np.sqrt(np.pi)*weights*x*x
    return np.sum(occupation*scale,axis=-1),np.sum(occupation*(1-occupation)*scale,axis=-1)


def hamiltonian(z,potential_eV,mass,states):
    """Dirichlet wavefunctions, centred finite differences, normalized in 1/m."""
    z=np.asarray(z); v=np.asarray(potential_eV)
    if z.ndim!=1 or len(z)<5 or v.shape!=z.shape or not np.all(np.isfinite(v)):
        raise ValueError("Hamiltonian needs a finite 1D potential on >=5 nodes")
    h=np.diff(z)
    if np.min(h)<=0 or not np.allclose(h,h[0],rtol=1e-9,atol=0) or mass<=0:
        raise ValueError("Use a uniform increasing grid and positive effective mass")
    if not isinstance(states,int) or not 1<=states<=len(z)-2:
        raise ValueError("Invalid number of eigenstates")
    kinetic=HBAR**2/(2*mass*M0*Q*h[0]**2)
    diagonal=2*kinetic+v[1:-1]
    off=np.full(len(z)-3,-kinetic)
    energy,interior=eigh_tridiagonal(diagonal,off,select="i",select_range=(0,states-1))
    wave=np.zeros((len(z),states)); wave[1:-1]=interior/np.sqrt(h[0])
    return energy,wave


class MOS:
    def __init__(self,device=Device(),spacing_nm=.05,states=32):
        if not np.isfinite(spacing_nm) or spacing_nm<=0:
            raise ValueError("Mesh spacing must be positive and finite")
        intervals=max(40,2*int(np.ceil(device.silicon_nm/spacing_nm/2)))
        if intervals>2400:
            raise ValueError("Teaching solver is limited to 2400 intervals")
        if not isinstance(states,int) or not 2<=states<intervals:
            raise ValueError("Use 2..intervals-1 eigenstates per valley ladder")
        self.device=device; self.states=states
        self.z=np.linspace(0,device.silicon_nm*1e-9,intervals+1)
        self.h=self.z[1]; self.weights=np.full(intervals+1,self.h)
        self.weights[[0,-1]]=self.h/2
        # A is scaled by eps/h. Its unknown is psi in volts.
        self.diag=np.full(intervals+1,2.0)
        self.diag[[0,-1]]=1+device.cox*self.h/EPS
        self.off=np.full(intervals,-1.0)
        self.charge_scale=Q*self.weights*self.h/EPS
        self.gate_scale=device.cox*self.h/EPS

    def apply(self,psi):
        out=self.diag*psi
        out[:-1]+=self.off*psi[1:]; out[1:]+=self.off*psi[:-1]
        return out

    def band_solve(self,rhs,extra_diagonal=0):
        band=np.zeros((3,len(self.z)))
        band[1]=self.diag+extra_diagonal
        band[0,1:]=self.off; band[2,:-1]=self.off
        return solve_banded((1,1),band,rhs)

    def gate_rhs(self,vg):
        rhs=np.zeros(len(self.z)); rhs[[0,-1]]=self.gate_scale*(vg-self.device.gate_offset_V)
        return rhs

    def classical_density(self,psi):
        d=self.device; value,derivative=fermi_half((psi-d.ec_reference_eV)/d.vt)
        prefactor=0.
        for _,gv,mz,md in VALLEYS:
            prefactor+=2*gv*(M0*KB*d.temperature_K/(2*np.pi*HBAR**2))**1.5*md*np.sqrt(mz)
        return prefactor*value,prefactor*derivative/d.vt

    def quantum_density(self,psi):
        d=self.device; n=np.zeros(len(self.z)); ladders=[]
        for name,gv,mz,md in VALLEYS:
            energy,wave=hamiltonian(self.z,d.ec_reference_eV-psi,mz,self.states)
            prefactor=gv*md*M0*KB*d.temperature_K/(np.pi*HBAR**2)
            populations=prefactor*np.logaddexp(0,-energy/d.vt)
            n+=wave**2@populations
            ladders.append(dict(name=name,degeneracy=gv,mass_z=mz,mass_dos=md,
                                energies_eV=energy,wavefunctions=wave,sheet_populations_m2=populations,
                                dNs_dEF_eV=prefactor/d.vt*expit(-energy/d.vt)))
        total=sum(np.sum(x["sheet_populations_m2"]) for x in ladders)
        tail=sum(np.sum(x["sheet_populations_m2"][-3:]) for x in ladders)/max(total,1.0)
        return n,ladders,float(tail)

    def residual(self,psi,n,vg):
        return self.apply(psi)+self.charge_scale*(self.device.na+n)-self.gate_rhs(vg)

    def solve(self,vg,model="quantum",initial=None,mixing=.6,tolerance=2e-11,max_outer=400):
        if not np.isfinite(vg) or not .05<=vg<=1.1:
            raise ValueError("Electron-only teaching branch: Vg=.05..1.1 V")
        if model not in ("quantum","classical") or not 0<mixing<=1 or tolerance<=0:
            raise ValueError("Invalid model, mixing or tolerance")
        psi=self.band_solve(self.gate_rhs(vg)-self.charge_scale*self.device.na) if initial is None else np.array(initial,dtype=float).copy()
        if psi.shape!=self.z.shape or not np.all(np.isfinite(psi)):
            raise ValueError("Initial potential must match grid and be finite")
        history=[]; d=self.device
        for iteration in range(max_outer):
            if model=="quantum":
                n,ladders,tail=self.quantum_density(psi)
            else:
                n,_=self.classical_density(psi); ladders=[]; tail=0.
            # Predictor is a preconditioner only; accepted residual uses recomputed
            # quantum density, not this local exponential response.
            trial=psi.copy()
            for inner in range(60):
                if model=="quantum":
                    exponent=(trial-psi)/d.vt
                    if np.max(abs(exponent))>80: raise RuntimeError("Predictor overflow")
                    pred=n*np.exp(exponent); derivative=pred/d.vt
                else: pred,derivative=self.classical_density(trial)
                residual=self.residual(trial,pred,vg)
                if np.max(abs(residual))<8e-16: break
                step=self.band_solve(-residual,self.charge_scale*derivative)
                factor=min(1.,.08/max(np.max(abs(step)),1e-30))
                baseline=np.linalg.norm(residual)
                for _ in range(25):
                    candidate=trial+factor*step
                    if model=="quantum":
                        candidate_n=n*np.exp((candidate-psi)/d.vt)
                    else: candidate_n,_=self.classical_density(candidate)
                    if np.linalg.norm(self.residual(candidate,candidate_n,vg))<baseline:
                        trial=candidate; break
                    factor/=2
                else: raise RuntimeError("Poisson predictor line search failed")
            else: raise RuntimeError("Poisson predictor did not converge")
            update=float(np.max(abs(trial-psi)))
            # Equivalent correction from the true charge, in V: guards against a
            # tiny mixed update masquerading as physical self-consistency.
            correction=float(np.max(abs(self.band_solve(self.residual(psi,n,vg)))))
            history.append(dict(iteration=iteration,predictor_update_V=update,true_poisson_correction_V=correction))
            if update<tolerance and correction<2*tolerance:
                if tail>1e-8: raise RuntimeError("Occupied-state tail too large; increase eigenstate count")
                return self.result(vg,model,psi,n,ladders,tail,history)
            psi+=mixing*(trial-psi)
        raise RuntimeError(f"Self-consistency failed at Vg={vg}, model={model}; correction={correction}")

    def result(self,vg,model,psi,n,ladders,tail,history):
        d=self.device; ns=float(self.weights@n)
        distance=np.minimum(self.z,self.z[-1]-self.z)
        centroid=float(self.weights@(n*distance)/max(ns,1e-100)*1e9)
        gate=d.cox*(vg-d.gate_offset_V-psi[[0,-1]])
        silicon=-Q*(self.weights@(n+d.na))
        denom=max(abs(silicon),1e-20)
        return dict(model=model,vg=vg,z=self.z.copy(),potential_V=psi.copy(),density_m3=n.copy(),
                    Ns_m2=ns,centroid_nm=centroid,gate_charge_C_m2=float(gate.sum()),
                    silicon_charge_C_m2=float(silicon),relative_gauss_error=float(abs(gate.sum()+silicon)/denom),
                    ladders=ladders,tail_fraction=tail,history=history,
                    residual_correction_V=history[-1]["true_poisson_correction_V"])

    def sweep(self,biases,model="quantum",mixing=.6):
        result=[]; previous=None
        for vg in biases:
            r=self.solve(float(vg),model,previous,mixing=mixing)
            result.append(r); previous=r["potential_V"]
        return result

    def capacitance(self,vg,model="quantum",step=.001,initial=None):
        if not np.isfinite(step) or step<=0: raise ValueError("Voltage step must be positive")
        a=self.solve(vg-step,model,initial); b=self.solve(vg+step,model,initial)
        return (b["gate_charge_C_m2"]-a["gate_charge_C_m2"])/(2*step)
