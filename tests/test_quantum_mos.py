"""Physics, statistics, independent references and convergence checks."""
from pathlib import Path
import sys,unittest
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from quantum_mos import MOS,Device,hamiltonian,fermi_half,Q,KB,HBAR,M0
from reference import box_levels,triangular_levels,finite_element_spectrum,coupled_root


class QuantumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=MOS();cls.q=cls.m.solve(.8);cls.c=cls.m.solve(.8,"classical")

    def test_box_energies_and_second_order(self):
        errors=[]
        for h in [.2,.1,.05]:
            z=np.linspace(0,10e-9,int(10/h)+1)
            e,_=hamiltonian(z,np.zeros(len(z)),.916,3)
            errors.append(np.max(abs(e/box_levels(10e-9,.916,3)-1)))
        self.assertGreater(np.log2(errors[-2]/errors[-1]),1.98)

    def test_triangular_airy_levels(self):
        z=np.linspace(0,20e-9,801)
        e,_=hamiltonian(z,z*1e9*.05,.916,3)
        self.assertLess(np.max(abs(e/triangular_levels(.05,.916,3)-1)),7e-5)

    def test_wavefunction_normalization_orthogonality_and_walls(self):
        for ladder in self.q["ladders"]:
            w=ladder["wavefunctions"]
            np.testing.assert_allclose(w.T@(self.m.weights[:,None]*w),np.eye(self.m.states),atol=2e-12)
            np.testing.assert_array_equal(w[[0,-1]],0)

    def test_sheet_density_matches_subband_population(self):
        total=sum(np.sum(a["sheet_populations_m2"]) for a in self.q["ladders"])
        self.assertLess(abs(total/self.q["Ns_m2"]-1),1e-12)

    def test_spin_valley_occupation_against_energy_integral(self):
        d=self.m.device
        for a in self.q["ladders"]:
            eta=-a["energies_eV"][0]/d.vt
            integral=quad(lambda t:expit(eta-t),0,max(eta,0)+60,epsabs=1e-12)[0]
            density=a["degeneracy"]*a["mass_dos"]*M0*KB*d.temperature_K/(np.pi*HBAR**2)*integral
            self.assertLess(abs(density/a["sheet_populations_m2"][0]-1),1e-10)

    def test_classical_fermi_integral_and_derivative(self):
        for eta in [-20,0,10,30,80]:
            f,df=fermi_half(np.array([eta]))
            end=np.sqrt(max(eta,0)+70)
            expected=4/np.sqrt(np.pi)*quad(lambda x:x*x*expit(eta-x*x),0,end,epsabs=1e-13,epsrel=1e-12)[0]
            deriv=4/np.sqrt(np.pi)*quad(lambda x:x*x*expit(eta-x*x)*(1-expit(eta-x*x)),0,end,epsabs=1e-13,epsrel=1e-12)[0]
            self.assertLess(abs(f[0]/expected-1),1e-7)
            self.assertLess(abs(df[0]/deriv-1),2e-7)

    def test_charge_balance_and_true_poisson_residual(self):
        for r in [self.q,self.c]:
            self.assertLess(r["relative_gauss_error"],1e-7)
            self.assertLess(r["residual_correction_V"],4e-11)

    def test_double_gate_symmetry_and_density_positive(self):
        for r in [self.q,self.c]:
            np.testing.assert_allclose(r["potential_V"],r["potential_V"][::-1],atol=1e-11)
            np.testing.assert_allclose(r["density_m3"],r["density_m3"][::-1],rtol=1e-9,atol=1e6)
            self.assertGreaterEqual(r["density_m3"].min(),0)

    def test_confinement_charge_and_centroid(self):
        self.assertLess(self.q["Ns_m2"],self.c["Ns_m2"])
        self.assertGreater(self.q["centroid_nm"],self.c["centroid_nm"])

    def test_independent_fem_spectrum(self):
        e,_=finite_element_spectrum(self.m.z,self.m.device.ec_reference_eV-self.q["potential_V"],.916)
        self.assertLess(abs(e[0]-self.q["ladders"][0]["energies_eV"][0]),1.3e-5)

    def test_alternate_coupled_root_algorithm(self):
        m=MOS(spacing_nm=.25);a=m.solve(.8);b=coupled_root(m,.8)
        self.assertLess(np.max(abs(a["potential_V"]-b["potential_V"])),1e-8)
        self.assertLess(abs(a["Ns_m2"]/b["Ns_m2"]-1),1e-7)

    def test_mixing_and_bias_path_independence(self):
        other=self.m.solve(.8,mixing=.35)
        reverse=self.m.sweep([1.,.9,.8])[-1]
        for r in [other,reverse]: self.assertLess(abs(r["Ns_m2"]/self.q["Ns_m2"]-1),2e-8)

    def test_states_and_mesh_convergence(self):
        more=MOS(states=48);n,_,_=more.quantum_density(self.q["potential_V"])
        self.assertLess(abs((more.weights@n)/self.q["Ns_m2"]-1),1e-10)
        fine=MOS(spacing_nm=.025).solve(.8)
        self.assertLess(abs(fine["Ns_m2"]/self.q["Ns_m2"]-1),2e-4)

    def test_gate_capacitance_and_voltage_step(self):
        a=self.m.capacitance(.8,initial=self.q["potential_V"])
        b=self.m.capacitance(.8,step=.0005,initial=self.q["potential_V"])
        self.assertGreater(a,0);self.assertLess(a,2*self.m.device.cox)
        self.assertLess(abs(a/b-1),2e-5)

    def test_occupation_tail_and_nonconvergence_are_rejected(self):
        with self.assertRaises(RuntimeError): MOS(states=4).solve(.8)
        with self.assertRaises(RuntimeError): self.m.solve(.8,max_outer=1)

    def test_input_guards(self):
        for kw in [{"silicon_nm":0},{"temperature_K":np.nan},{"oxide_nm":-1}]:
            with self.assertRaises(ValueError): Device(**kw)
        with self.assertRaises(ValueError): MOS(spacing_nm=0)
        with self.assertRaises(ValueError): self.m.solve(.8,mixing=0)
        with self.assertRaises(ValueError): self.m.solve(-1)
        with self.assertRaises(ValueError): self.m.capacitance(.8,step=0)


if __name__=="__main__":unittest.main()
