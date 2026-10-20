"""Generate eight figures with raw numerical data and verification metrics."""
from pathlib import Path
from dataclasses import replace,asdict
import csv,json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from quantum_mos import MOS,Device,VALLEYS,hamiltonian,Q
from reference import finite_element_spectrum,box_levels,triangular_levels,coupled_root

plt.rcParams.update({"font.size":10,"axes.grid":True,"grid.alpha":.22,"figure.dpi":120,"savefig.dpi":170})


def csv_out(path,rows):
    rows=list(rows)
    with Path(path).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def save(fig,out,name):
    fig.savefig(out/(name+".png"),bbox_inches="tight");plt.close(fig)


def row(r):
    populations=[float(np.sum(v["sheet_populations_m2"])) for v in r["ladders"]]
    return dict(model=r["model"],Vg_V=r["vg"],Ns_cm2=r["Ns_m2"]/1e4,
                nearest_interface_centroid_nm=r["centroid_nm"],
                gate_charge_C_m2=r["gate_charge_C_m2"],
                residual_correction_V=r["residual_correction_V"],
                relative_gauss_error=r["relative_gauss_error"],iterations=len(r["history"]),
                Delta2_fraction=populations[0]/r["Ns_m2"] if populations else None,
                Delta2_ground_eV=r["ladders"][0]["energies_eV"][0] if populations else None,
                Delta4_ground_eV=r["ladders"][1]["energies_eV"][0] if populations else None)


def baseline(out,spacing=.05):
    m=MOS(spacing_nm=spacing); biases=np.linspace(.1,1.0,37)
    quantum=m.sweep(biases);classical=m.sweep(biases,"classical")
    data=[]
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for results,label in [(quantum,"Self-consistent quantum"),(classical,"Classical 3D Fermi-Dirac")]:
        charge=np.array([r["gate_charge_C_m2"] for r in results])
        capacitance=np.gradient(charge,biases,edge_order=2)
        axs[0].semilogy(biases,[r["Ns_m2"]/1e4 for r in results],label=label)
        axs[1].plot(biases[1:-1],capacitance[1:-1]/(2*m.device.cox),label=label)
        for r,c in zip(results,capacitance):
            value=row(r);value["dQg_dVg_F_m2"]=float(c);data.append(value)
    axs[0].axhline(1e12,color="gray",ls="--",label="Sheet-charge criterion")
    axs[0].set(xlabel="Tied gate voltage (V)",ylabel="Electron sheet density (cm$^{-2}$)",title="Same bands, reservoir and Fermi-Dirac statistics")
    axs[1].axhline(1,color="gray",ls="--",label="Two oxides in parallel")
    axs[1].set(xlabel="Tied gate voltage (V)",ylabel="d(Qg,top+Qg,bottom)/dVg / (2 Cox)",title="Equilibrium quasi-static gate capacitance")
    for ax in axs:ax.legend(fontsize=8)
    save(fig,out,"charge_capacitance");csv_out(out/"bias_sweep.csv",data)

    profile_rows=[];selected={}
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,vg in zip(axs,[.4,.8]):
        qr=m.solve(vg);cr=m.solve(vg,"classical");selected[vg]=(qr,cr)
        for r,label in [(qr,"Quantum"),(cr,"Classical FD")]:
            ax.plot(m.z*1e9,r["density_m3"]/1e6,label=label)
            profile_rows.extend(dict(Vg_V=vg,model=r["model"],z_nm=z*1e9,potential_V=p,n_cm3=n/1e6)
                                for z,p,n in zip(m.z,r["potential_V"],r["density_m3"]))
        ax.set(xlabel="Silicon z (nm)",ylabel="Electron density (cm$^{-3}$)",title=f"Vg={vg:.2f} V; infinite oxide barriers")
        ax.legend()
    save(fig,out,"density_profiles");csv_out(out/"density_profiles.csv",profile_rows)

    qr,cr=selected[.8];np.savez_compressed(out/"quantum_state.npz",z_m=m.z,
          potential_V=qr["potential_V"],density_m3=qr["density_m3"],
          Delta2_energies_eV=qr["ladders"][0]["energies_eV"],
          Delta4_energies_eV=qr["ladders"][1]["energies_eV"],
          Delta2_wavefunctions_per_sqrt_m=qr["ladders"][0]["wavefunctions"],
          Delta4_wavefunctions_per_sqrt_m=qr["ladders"][1]["wavefunctions"])
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    spectra=[]
    for ax,ladder in zip(axs,qr["ladders"]):
        ax.plot(m.z*1e9,m.device.ec_reference_eV-qr["potential_V"],color="black",lw=2,label="Conduction-band edge")
        for i in range(3):
            energy=ladder["energies_eV"][i];w=ladder["wavefunctions"][:,i]**2
            ax.axhline(energy,color=f"C{i}",ls="--",alpha=.6)
            ax.plot(m.z*1e9,energy+.025*w/w.max(),color=f"C{i}",label=f"E{i+1} + scaled probability")
        ax.axhline(0,color="gray",ls=":",label="Reservoir EF=0")
        ax.set(xlabel="z (nm)",ylabel="Energy (eV)",title=f"{ladder['name']}: mz={ladder['mass_z']:.3f} m0")
        ax.legend(fontsize=7)
        spectra.extend(dict(valley=ladder["name"],state=i+1,energy_eV=float(e),Ns_cm2=float(ns/1e4))
                       for i,(e,ns) in enumerate(zip(ladder["energies_eV"],ladder["sheet_populations_m2"])))
    save(fig,out,"bands_wavefunctions");csv_out(out/"subbands.csv",spectra)

    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for r,label in [(quantum,"Quantum"),(classical,"Classical FD")]:
        axs[0].plot(biases,[a["centroid_nm"] for a in r],label=label)
    axs[0].set(xlabel="Tied gate voltage (V)",ylabel="Mean distance to nearest interface (nm)",title="Global mean z alone would always be tSi/2")
    axs[0].legend()
    fraction=[np.sum(a["ladders"][0]["sheet_populations_m2"])/a["Ns_m2"] for a in quantum]
    axs[1].plot(biases,fraction,label="Delta2 population fraction")
    axs[1].axhline(2/6,color="gray",ls="--",label="Valley count alone: 2/6")
    axs[1].set(xlabel="Tied gate voltage (V)",ylabel="Fraction of total electron sheet density",title="Confinement and in-plane DOS compete",ylim=(0,1))
    axs[1].legend(fontsize=8)
    save(fig,out,"centroid_valleys")
    return m,qr,cr,quantum,classical


def criterion(model,kind,target=1e12):
    cache=None
    def f(v):
        nonlocal cache
        r=model.solve(v,kind,cache);cache=r["potential_V"]
        return np.log(r["Ns_m2"]/(target*1e4))
    return float(brentq(f,.05,1.1,xtol=2e-7))


def parameters(out,spacing=.05):
    rows=[]
    for tsi in [3,4,6,8,10,16]:
        m=MOS(Device(silicon_nm=tsi),spacing_nm=spacing)
        qc=criterion(m,"quantum");cc=criterion(m,"classical")
        rows.append(dict(silicon_nm=tsi,oxide_nm=2.,temperature_K=300.,
                         quantum_criterion_V=qc,classical_criterion_V=cc,shift_mV=1000*(qc-cc)))
    oxide_rows=[]
    for tox in [1,2,4]:
        m=MOS(Device(oxide_nm=tox),spacing_nm=spacing)
        qc=criterion(m,"quantum");cc=criterion(m,"classical")
        oxide_rows.append(dict(silicon_nm=10.,oxide_nm=tox,quantum_criterion_V=qc,
                               classical_criterion_V=cc,shift_mV=1000*(qc-cc)))
    csv_out(out/"thickness_criterion.csv",rows);csv_out(out/"oxide_criterion.csv",oxide_rows)
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    axs[0].plot([r["silicon_nm"] for r in rows],[r["shift_mV"] for r in rows],"o-")
    axs[0].set(xlabel="Silicon thickness (nm)",ylabel="Quantum - classical criterion voltage (mV)",title="Fixed Ns=1e12 cm$^{-2}$; tox=2 nm")
    axs[1].plot([r["oxide_nm"] for r in oxide_rows],[r["quantum_criterion_V"] for r in oxide_rows],"o-",label="Quantum")
    axs[1].plot([r["oxide_nm"] for r in oxide_rows],[r["classical_criterion_V"] for r in oxide_rows],"o-",label="Classical FD")
    axs[1].set(xlabel="Each oxide thickness (nm)",ylabel="Sheet-charge criterion gate voltage (V)",title="Same criterion and reservoir alignment")
    axs[1].legend()
    save(fig,out,"geometry_criterion")
    return rows,oxide_rows


def verification(out,spacing=.05):
    mesh=[];spectral=[]
    for h in [2*spacing,spacing,spacing/2]:
        m=MOS(spacing_nm=h);r=m.solve(.8)
        cap=m.capacitance(.8,initial=r["potential_V"])
        fem,_=finite_element_spectrum(m.z,m.device.ec_reference_eV-r["potential_V"],.916)
        fd=r["ladders"][0]["energies_eV"][:6]
        mesh.append(dict(spacing_nm=m.h*1e9,Ns_cm2=r["Ns_m2"]/1e4,centroid_nm=r["centroid_nm"],
                         ground_energy_eV=fd[0],capacitance_F_m2=cap,
                         FD_FEM_ground_difference_meV=1000*abs(fd[0]-fem[0])))
        z=np.linspace(0,20e-9,int(round(20/h))+1)
        box,_=hamiltonian(z,np.zeros(len(z)),.916,3)
        triangle,_=hamiltonian(z,.05*z*1e9,.916,3)
        spectral.append(dict(spacing_nm=h,box_max_relative_error=float(np.max(abs(box/box_levels(20e-9,.916,3)-1))),
                             airy_max_relative_error=float(np.max(abs(triangle/triangular_levels(.05,.916,3)-1)))))
    csv_out(out/"mesh_convergence.csv",mesh);csv_out(out/"exact_spectra.csv",spectral)

    base=MOS(spacing_nm=spacing);r=base.solve(.8)
    counts=[]
    for states in [4,8,16,32,48]:
        model=MOS(spacing_nm=spacing,states=states)
        n,_,tail=model.quantum_density(r["potential_V"])
        ns=model.weights@n
        counts.append(dict(states_per_ladder=states,Ns_frozen_potential_cm2=ns/1e4,
                           relative_difference_vs_32=float(ns/r["Ns_m2"]-1),retained_last_three_fraction=tail))
    steps=[]
    for dv in [.002,.001,.0005]:
        steps.append(dict(dV_V=dv,gate_capacitance_F_m2=base.capacitance(.8,step=dv,initial=r["potential_V"])))
    csv_out(out/"state_truncation.csv",counts);csv_out(out/"capacitance_steps.csv",steps)
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
    hs=[a["spacing_nm"] for a in mesh]
    axs[0].plot(hs,100*(np.array([a["Ns_cm2"] for a in mesh])/mesh[-1]["Ns_cm2"]-1),"o-")
    axs[0].set(xlabel="Mesh spacing (nm)",ylabel="Sheet density relative to finest (%)",
               title=f"Mesh comparison; finest={hs[-1]:g} nm")
    axs[1].semilogy([a["states_per_ladder"] for a in counts],
        [max(abs(a["relative_difference_vs_32"]),1e-15) for a in counts],"o-")
    axs[1].set(xlabel="Eigenstates per valley ladder",ylabel="Relative sheet-charge difference",title="Frozen potential; 32-state reference")
    axs[2].plot(np.array([a["dV_V"] for a in steps])*1000,[a["gate_capacitance_F_m2"] for a in steps],"o-")
    axs[2].ticklabel_format(axis="y",style="plain",useOffset=False)
    axs[2].set(xlabel="Central-difference step (mV)",ylabel="Tied-gate capacitance (F/m$^2$)",title="Every bias re-converges quantum charge")
    save(fig,out,"convergence")
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for key,label in [("box_max_relative_error","Infinite box, first 3 levels"),("airy_max_relative_error","Triangular well, first 3 levels")]:
        axs[0].loglog(hs,[a[key] for a in spectral],"o-",label=label)
    axs[0].set_xticks(hs,labels=[f"{v:g}" for v in hs]);axs[0].tick_params(axis="x",which="minor",labelbottom=False,bottom=False)
    axs[0].set(xlabel="Grid spacing (nm)",ylabel="Maximum relative energy error",title="Two analytically known potentials")
    axs[0].legend(fontsize=8)
    axs[1].loglog(hs,[a["FD_FEM_ground_difference_meV"] for a in mesh],"o-")
    axs[1].set_xticks(hs,labels=[f"{v:g}" for v in hs]);axs[1].tick_params(axis="x",which="minor",labelbottom=False,bottom=False)
    axs[1].set(xlabel="Grid spacing (nm)",ylabel="Ground-state difference (meV)",title="Frozen SP potential: FD vs independent P1 FEM")
    save(fig,out,"spectrum_validation")
    return dict(mesh=mesh,exact_spectra=spectral,capacitance_steps=steps,
                state_truncation=counts,exact_box_final_order=float(np.log(spectral[-2]["box_max_relative_error"]/spectral[-1]["box_max_relative_error"])/np.log(2)))


def feedback(out,spacing=.05):
    m=MOS(spacing_nm=spacing); classical=m.solve(.8,"classical");full=m.solve(.8)
    one_n,_,_=m.quantum_density(classical["potential_V"])
    one_ns=float(m.weights@one_n)
    coarse=MOS(spacing_nm=.25);primary=coarse.solve(.8);alternate=coupled_root(coarse,.8)
    mixrows=[];traces={}
    for mix in [.35,.6,.9]:
        r=m.solve(.8,mixing=mix);traces[str(mix)]=r["history"]
        mixrows.append(dict(mixing=mix,iterations=len(r["history"]),Ns_cm2=r["Ns_m2"]/1e4,
                            max_potential_difference_V=float(np.max(abs(r["potential_V"]-full["potential_V"])))))
    csv_out(out/"mixing_check.csv",mixrows)
    (out/"iteration_traces.json").write_text(json.dumps(traces,indent=2)+"\n")
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for key,trace in traces.items():
        axs[0].semilogy([a["iteration"] for a in trace],[max(a["true_poisson_correction_V"],1e-15) for a in trace],label=f"Mix={key}")
    axs[0].set(xlabel="Outer iteration",ylabel="True-charge Poisson correction (V)",title="Check residual before mixing")
    axs[0].legend()
    axs[1].plot(m.z*1e9,one_n/1e6,label="One-shot quantum on classical potential")
    axs[1].plot(m.z*1e9,full["density_m3"]/1e6,label="Self-consistent quantum")
    axs[1].set(xlabel="Silicon z (nm)",ylabel="Electron density (cm$^{-3}$)",title="Quantum charge must feed back into Poisson")
    axs[1].legend(fontsize=8)
    save(fig,out,"self_consistency")
    return dict(one_shot_Ns_cm2=one_ns/1e4,self_consistent_Ns_cm2=full["Ns_m2"]/1e4,
                one_shot_relative_difference=one_ns/full["Ns_m2"]-1,
                alternate_root_grid_nm=.25,alternate_root_evaluations=alternate["evaluations"],
                alternate_root_potential_difference_V=float(np.max(abs(primary["potential_V"]-alternate["potential_V"]))),
                alternate_root_relative_Ns_difference=primary["Ns_m2"]/alternate["Ns_m2"]-1,mixing=mixrows)


def run(output,spacing=.05):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    model,qr,cr,quantum,classical=baseline(out,spacing)
    geometry,oxide=parameters(out,spacing)
    validation=verification(out,spacing);validation["feedback"]=feedback(out,spacing)
    validation["maximum_sweep_gauss_error"]=max(r["relative_gauss_error"] for r in quantum+classical)
    validation["maximum_sweep_poisson_correction_V"]=max(r["residual_correction_V"] for r in quantum+classical)
    capq=model.capacitance(.8,initial=qr["potential_V"])
    capc=model.capacitance(.8,"classical",initial=cr["potential_V"])
    basecriterion=next(r for r in geometry if r["silicon_nm"]==10)
    summary=dict(device=asdict(model.device),spacing_nm=model.h*1e9,states_per_ladder=model.states,
                 at_0p8V_quantum=row(qr),at_0p8V_classical=row(cr),criterion_1e12=basecriterion,
                 gate_capacitance_quantum_F_m2=capq,gate_capacitance_classical_F_m2=capc,
                 two_oxide_limit_F_m2=2*model.device.cox,
                 frozen_spectrum_quantum_capacitance_F_m2=float(Q*sum(np.sum(a["dNs_dEF_eV"]) for a in qr["ladders"])),
                 scope="Equilibrium electron-only infinite-barrier thin-Si Hartree model; no holes, transport, exchange-correlation or experimental data.")
    for name,data in [("summary.json",summary),("validation.json",validation)]:
        (out/name).write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2));print(f"Saved eight figures and reproducible data to {out}")
    return summary,validation
