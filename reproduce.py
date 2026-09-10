"""Regenerate every revised numerical result and figure from explicit inputs."""
import csv
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from diagnostics import (cube_integrand, cube_increment_direct, mean_se,
                         energy_coefficient, mode_count, universal_bound,
                         abc_flow, abc_velocity)

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'
SEED=20260910
N=4096


def write_csv(name,rows):
    with (RESULTS/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)


def main(output_dir=ROOT):
    global RESULTS, FIGURES
    output_dir = Path(output_dir).resolve()
    RESULTS = output_dir/'results'; FIGURES = output_dir/'figures'
    RESULTS.mkdir(parents=True,exist_ok=True);FIGURES.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(SEED)
    x=rng.random((N,3));np.save(RESULTS/'sample_points.npy',x)
    rows=[]
    for name,g,exact in [
        ('translation',np.full(N,np.log(3-2*np.cos(2*np.pi*.123))),np.log(3-2*np.cos(2*np.pi*.123))),
        ('cat_map',np.log(3-2*np.cos(2*np.pi*(x[:,0]+x[:,1]))),np.log((3+np.sqrt(5))/2))]:
        est,se=mean_se(g)
        rows.append(dict(case=name,exact=float(exact),estimate=est,se=se,N=N,seed=SEED))
        np.save(RESULTS/(name+'_integrand.npy'),g)
    write_csv('benchmarks.csv',rows)

    limits=[]
    for k in [1,2,4,8,16,32,64,128]:
        d=np.tile([.123,0,0],(1,1))
        s=float(cube_integrand(d,k)[0]);normalized=float(cube_integrand(d,k,True)[0])
        limits.append(dict(K=k,m=mode_count(k),S=s,S_over_K2=s/k**2,
                           bound_over_K2=universal_bound(k)/k**2,
                           S_over_log_m=s/np.log(mode_count(k)),normalized=normalized))
    write_csv('cutoff_translation.csv',limits)

    violations=[]
    for time,vals in [(0,[.12,.11,.11]),(2.5,[.84,.92,1.03]),(4.5,[3.21,4.87,7.14])]:
        for k,value in zip([2,4,8],vals):
            bound=universal_bound(k)/k**2
            violations.append(dict(time=time,K=k,original_claim=value,upper_bound=bound,
                                    incompatible=bool(value>bound)))
    write_csv('original_table_bound_audit.csv',violations)

    small=[];k=4
    sample_limit=4*np.pi**2*energy_coefficient(k)*np.sum(abc_velocity(x)**2,axis=1)
    exact_limit=4*np.pi**2*energy_coefficient(k)*3
    for lag in [2e-3,1e-3,5e-4,2e-4,1e-4,5e-5,2e-5]:
        y=abc_flow(x,lag,steps=16)
        g=cube_integrand(y-x,k)/lag**2
        value,se=mean_se(g)
        small.append(dict(h=lag,estimate=value,se=se,exact_limit=exact_limit,
                          sampled_limit=float(sample_limit.mean()),N=N))
    write_csv('small_time_abc.csv',small)

    abc=[];last=None
    for steps in [64,128,256]:
        y=abc_flow(x,1,steps=steps)
        g=cube_integrand(y-x,4);s,se=mean_se(g)
        difference=0.0 if last is None else float(np.max(np.linalg.norm(y-last,axis=1)))
        abc.append(dict(steps=steps,lag=1,viscosity=.01,start_time=0,K=4,
                        S=s,se=se,max_position_change_from_previous=difference,N=N))
        last=y
    np.save(RESULTS/'abc_final_positions.npy',last)
    np.save(RESULTS/'abc_final_integrand.npy',g)
    write_csv('abc_step_refinement.csv',abc)

    y=(np.arange(2**15)+.5)/2**15
    shears=[]
    for n in [1,3,9,27,81]:
        delta=np.zeros((len(y),3));delta[:,0]=.2*np.sin(2*np.pi*n*y)
        s=float(cube_integrand(delta,4).mean())
        gradient=float(np.mean((.2*2*np.pi*n*np.cos(2*np.pi*n*y))**2))
        shears.append(dict(n=n,h=.2,K=4,S=s,vorticity_sup=2*np.pi*n,
                           deformation_L2_squared=gradient,
                           deformation_exact=.2**2*(2*np.pi*n)**2/2))
    write_csv('shear_counterexample.csv',shears)

    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.bbox':'tight','pdf.fonttype':42,'ps.fonttype':42})
    fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    ks=np.array([v['K'] for v in limits])
    ax[0].loglog(ks,[v['S_over_K2'] for v in limits],'o-',label='Translation')
    ax[0].loglog(ks,[v['bound_over_K2'] for v in limits],'--',label='Universal upper bound')
    ax[0].set(xlabel='Fourier cutoff K',ylabel=r'$S_K/K^2$',title='(a) Universal decay');ax[0].legend()
    ax[1].semilogx(ks,[v['normalized'] for v in limits],'o-',label='Normalized inside log')
    ax[1].axhline(np.log(3),ls='--',color='k',label=r'$\log 3$')
    ax[1].set(xlabel='Fourier cutoff K',ylabel=r'$\overline{S}_K$',title='(b) Fixed-map saturation');ax[1].legend()
    fig.savefig(FIGURES/'cutoff_limits.pdf');plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.6,3.7),layout='constrained')
    hs=np.array([r['h'] for r in small]);vals=np.array([r['estimate'] for r in small])/exact_limit
    ax.errorbar(hs,vals,yerr=2*np.array([r['se'] for r in small])/exact_limit,fmt='o-',capsize=3,
                label='MC estimate with 2 SE')
    ax.axhline(1,color='k',ls='--',label='Exact spatial limit')
    ax.axhline(sample_limit.mean()/exact_limit,color='tab:orange',ls=':',label='Same-sample limit')
    ax.set(xscale='log',xlabel='Time lag h',ylabel=r'$S_h/(h^2 C_K\|u(0)\|_2^2)$',
           title='Viscous ABC velocity: small-time energy limit');ax.legend(fontsize=8)
    fig.savefig(FIGURES/'small_time.pdf');plt.close(fig)

    fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    ns=[v['n'] for v in shears]
    ax[0].semilogx(ns,[v['S'] for v in shears],'o-');ax[0].set(xlabel='Shear frequency n',
       ylabel=r'$S_4(\Phi_h)$',title='(a) Value diagnostic is independent of n')
    ax[1].loglog(ns,[v['deformation_L2_squared'] for v in shears],'o-',label='Quadrature')
    ax[1].loglog(ns,[v['deformation_exact'] for v in shears],'--',label='Exact quadratic growth')
    ax[1].set(xlabel='Shear frequency n',ylabel=r'$\|D\Phi_h-I\|_{L^2}^2$',
              title='(b) Derivative probe detects deformation');ax[1].legend()
    fig.savefig(FIGURES/'shear_blindness.pdf');plt.close(fig)

    # TeX table written from the actual computed rows, with exact inputs preserved.
    tex=[r'\begin{tabular}{lrrr}',r'\toprule',r'Case & Exact value & Estimate & SE\\',r'\midrule']
    for row in rows:
        label={'translation':'Translation','cat_map':'Cat map'}[row['case']]
        tex.append(f"{label} & {row['exact']:.9f} & {row['estimate']:.9f} & {row['se']:.9f}"+r'\\')
    tex += [r'\bottomrule',r'\end{tabular}']
    (output_dir/'benchmark_table.tex').write_text('\n'.join(tex)+'\n')
    metadata={'seed':SEED,'N':N,'python':sys.version,'platform':platform.platform(),
              'numpy':np.__version__,'matplotlib':matplotlib.__version__,
              'scope':'Fresh audit calculations; not a reproduction of unavailable original code.',
              'sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in [ROOT/'diagnostics.py',ROOT/'reproduce.py',ROOT/'tests/test_diagnostics.py',
                                  ROOT/'verify_reproduction.py', ROOT/'requirements.txt']}}
    (RESULTS/'environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({'benchmarks':rows,'abc':abc,'small_time_last':small[-1],
                      'incompatible_original_cells':sum(v['incompatible'] for v in violations),
                      'shear_S_range':[min(v['S'] for v in shears),max(v['S'] for v in shears)]},indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT,
                        help='Destination for results, figures and the TeX table; defaults to the project root.')
    main(parser.parse_args().output_dir)
