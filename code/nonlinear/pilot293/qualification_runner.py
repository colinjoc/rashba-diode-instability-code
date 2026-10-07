"""One integrated finite-beta supporting qualification. Import executes no solve.
Frozen complete schedule, connected branch, source, tail, mixed and C1/C2 checks.
The first insufficient inequality is saved; no retries or parameter replacement.
"""
import importlib.util,json,math,resource,time,hashlib,sys,os,tempfile
from pathlib import Path
import numpy as np


def publish_checkpoint(path, raw, limit):
    """Preserve the last complete bounded checkpoint if publication fails."""
    if len(raw.encode()) > limit:
        raise RuntimeError('frozen retained-output cap')
    path = Path(path)
    fd, name = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            handle.write(raw); handle.flush(); os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name): os.unlink(name)

def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def hierarchy(values,noise):
    d1,d2=values[1]-values[0],values[2]-values[1];b0,b1,b2=noise;n=b0+2*b1+b2
    ok=abs(d2)<=.75*abs(d1)+n and (d1*d2>=0 or min(abs(d1),abs(d2))<=n)
    return dict(d1=d1,d2=d2,noise=n,pass_=bool(ok),error=4*(abs(d2)+b1+b2)+b2)

def run(package,output,guard):
    package,output=Path(package),Path(output);output.mkdir(exist_ok=False)
    p=json.loads((package/'QUALIFICATION_PROTOCOL.json').read_text())
    sol=load(package/'nonlinear_response.py','nonlinear_response');bd=load(package/'correction_bounds.py','finite_beta_bounds')
    anchorpath=package/'LINEAR271.json'
    if hashlib.sha256(anchorpath.read_bytes()).hexdigest()!=p['linear_anchor_sha256']:raise ValueError('linear anchor changed')
    anchor=json.loads(anchorpath.read_text());sol.require(anchor['qualification_pass'],'accepted linear anchor')
    f=p['frozen_state'];T,D,q,h=f['x'];rt,rD,rq,rh=p['source_box_radii'];K=f['K'][0];floor=p['arithmetic_floor']
    va=complex(f['v']['real'][0],f['v']['imag'][0]);vp=complex(f['v']['real'][1],f['v']['imag'][1])
    eta=np.array([(va-1j*vp).conjugate(),0.,va+1j*vp],dtype=np.complex128);E=sol.norm(eta)
    start,cpu=time.monotonic(),time.process_time();nodes=0;branchcells=0;validationnodes=0;rows={};context={};frequency_counts={};last_completed=None
    report=dict(schema='rashba-finite-beta-two-slab-control-cost-pilot-v1',scope=p['scope'],outcome='incomplete',qualification_pass=False,
        whole_qualification_pass=False,target_release=False,parent_answered=False,target_nodes=0,endpoint_confirmation_nodes=0,rows=[],failed_checks=[],
        uncomputed_frequencies=p['uncomputed_frequencies'],completed_frequencies=[],partial_sums={},data_role=p['data_role'],original_criterion=f,primitive_empirical=True,directed_rounding=False,whole_memory_certified=False,
        analytic_enclosures_conditional_on_empirical_binary64_allowances=True,
        exact_finite_Fourier_mean_without_collocation=True,
        accepted_linear_source_once=p['accepted_linear_source_bound'],nonlinear_source_allowance=p['nonlinear_source_allowance'],
        paired_odd_cancellation_before_absolute_error=True,matched_beta0='U0,V0,F0 computed within each new Taylor expansion; not replay271')
    def save():
        report['last_completed_node']=last_completed
        report['completed_nodes_per_frequency']=frequency_counts
        report['measurements']=dict(calculation_cpu_seconds=time.process_time()-cpu,calculation_wall_seconds=time.monotonic()-start,
            actual_two_sector_expansion_nodes=nodes,continuous_validation_expansion_nodes=validationnodes,
            connected_signed_beta_cells=branchcells,process_lifetime_maxrss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
            memory_scope='SELF lifetime maximum; native sampled PGID remains separate/incomplete; remote model memory excluded')
        raw=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n'
        if len(raw.encode())>p['artifact_bytes']:raise RuntimeError('frozen retained-output cap')
        publish_checkpoint(output/'QUALIFICATION_RESULT.json', raw, p['artifact_bytes'])
    def expansion(n,lam,theta,M,J,angle_radius=0.,validation=False):
        nonlocal nodes,validationnodes,context,branchcells,last_completed
        guard();u=0. if theta in (math.pi/2,3*math.pi/2) else math.cos(theta)
        # Exact global cos-cell enclosure, tight at exact grazing/extrema.
        ru=min(2.,abs(math.sin(theta))*angle_radius+.5*angle_radius*angle_radius)
        du=abs(q/2-lam*h)*ru+(rq/2+rh)*(abs(u)+ru)
        z=complex((2*n+1)*math.pi*T,(q/2-lam*h)*u);rz=(2*n+1)*math.pi*rt+du
        O,a,rO,ra=sol.background(z,D,rz,rD,floor);degree=max(M,2*J) if n>=p['analytic_tail_split'] else M
        context=dict(n=n,helicity=lam,theta=theta,profile_degree=degree,force_degree=2*J,angle_radius=angle_radius,source_box=p['source_box_radii'])
        report['current_attempt']=dict(context);save()
        U,eU=sol.sector(eta,a,O,D,u*K,degree,rO=rO,ra=ra,rD=rD,ruK=K*ru,floor=floor,guard=guard)
        V,eV=sol.sector(np.conj(eta[::-1]),a,O,D,-u*K,degree,rO=rO,ra=ra,rD=rD,ruK=K*ru,floor=floor,guard=guard)
        F,eF,coeff,identity,hpoly=sol.force_series(U,V,a,eta,2*J,eU=eU,eV=eV,ra=ra,floor=floor,guard=guard)
        # Formal coefficient normalization identity, all Fourier modes retained.
        defect=math.fsum(sol.norm(x)*p['epsilon'][0]**k for k,x in enumerate(identity[:len(F)]))
        sol.require(defect<=p['normalization_max'],'full Fourier coefficient normalized identity',defect=defect)
        residualA=sol.riccati_residual(U,eta,a,O,D,u*K,guard);residualB=sol.riccati_residual(V,np.conj(eta[::-1]),a,O,D,-u*K,guard)
        defectR=max(math.fsum(sol.norm(x)*p['epsilon'][0]**k for k,x in enumerate(rr[:len(U)])) for rr in (residualA,residualB))
        sol.require(defectR<=p['residual_max'],'both Taylor coefficient recurrence residuals',full_low_order_residual=defectR)
        force_omitted=max((len(x)//2 for x in identity),default=0)
        values=[]
        if n<p['analytic_tail_split']:
            rA=bd.residual_intervals(U,eU,eta,a,O,D,u*K,ra,rO,rD,K*ru,guard)
            rB=bd.residual_intervals(V,eV,np.conj(eta[::-1]),a,O,D,-u*K,ra,rO,rD,K*ru,guard)
            idI=bd.force_identity_intervals(U,V,F,eU,eV,eF,a,ra,guard)
            B=p['epsilon'][0];count=p['signed_beta_cells'];radius=B/count
            for j in range(count):
                b=-B+(2*j+1)*radius;context.update(beta_cell=[b-radius,b+radius])
                A=bd.branch_remainder(U,eU,eta,a,O,D,u*K,rA,b,radius,ra=ra,rO=rO,rD=rD,ruK=K*ru,H=p['preconditioner_modes'],cap=p['common_profile_remainder_cap'],floor=floor,guard=guard)
                BB=bd.branch_remainder(V,eV,np.conj(eta[::-1]),a,O,D,-u*K,rB,b,radius,ra=ra,rO=rO,rD=rD,ruK=K*ru,H=p['preconditioner_modes'],cap=p['common_profile_remainder_cap'],floor=floor,guard=guard)
                ff=bd.force_remainder(U,V,F,idI,eta,a,ra,b,radius,eU,eV,eF,A['radius'],BB['radius'],floor)
                branchcells+=1;values.append(ff['error'])
            remainder=max(values);primitive=floor*(1+sum(abs(c)*B**k for k,c in enumerate(coeff)))
        else:
            # The complex disk proof provides the branch and all omitted force
            # orders here. Exact finite coefficients through2J are required.
            cf=bd.tail_cf(T,rt,D,rD,p['analytic_beta_radius'],E,n,floor)
            B=p['epsilon'][0];ratio=B/p['analytic_beta_radius'];first=2*(len(F)//2)
            # An even-degree force polynomial's next even order is degree+2.
            first=2*((len(F)-1)//2+1)
            remainder=cf['Cf']/(cf['first_omitted']**3)*ratio**first/(1-ratio*ratio)+floor
            primitive=floor*(1+sum(abs(c)*B**k for k,c in enumerate(coeff)))
        src=[E*x for x in eF]
        nodes+=1;validationnodes+=int(validation)
        frequency_counts[str(n)]=frequency_counts.get(str(n),0)+1
        last_completed=dict(n=n,helicity=lam,theta=theta,angle_radius=angle_radius,validation=validation,normalization_defect=defect,low_order_residual=defectR)
        save() # after all required local checks; no uncompleted-node inflation
        return dict(coeff=coeff,source_coeff_error=src,remainder=remainder,primitive=primitive,
            low_order_residual=defectR,normalization_defect=defect,force_residual_max_mode=force_omitted,
            same_z=True,exact_grazing=u==0.)
    def moments(Nt):
        x,w=np.polynomial.legendre.leggauss(Nt);tt=(x+1)/2;ww=w/2
        sol.require(np.all(ww>0) and abs(float(np.sum(ww*tt))-.5)<=1e-13,'positive frozen Gauss path moment')
        return [float(np.sum(ww*tt**(k+1))) for k in range(2*max(p['levels'][3])+1)]
    try:
        sol.require(abs(K)>0 and f['K'][1]==0,'exact source-compatible fixed K')
        sol.require(abs(anchor['source_bridge']['new_conservative_bound']-p['accepted_linear_source_bound'])<=1e-20,'exact accepted linear source once')
        report['complete_paired_infinite_tail_envelopes']={str(512):{str(e):bd.paired_tail(T,rt,D,rD,e,E,512,p['analytic_beta_radius'],floor) for e in p['epsilon']}}
        report['tail_scope']='Complete conditional tail beyond512 only; missing finite slabs remain uncomputed.'
        save()
        M,Ntheta,J,Nt=p['slab_configuration']; moment=moments(Nt)
        gg=[float(t) for t in (np.polynomial.legendre.leggauss(Nt)[0]+1)/2]
        for n in p['frequencies']:
            guard();freq_start=time.monotonic();freq_cpu=time.process_time();before_nodes=nodes;before_validation=validationnodes
            sums={str(e):dict(correction=0.,source_force_error=0.,primitive_error=0.,residual_force_error=0.,path_polynomial_error=0.) for e in p['epsilon']}
            report['partial_sums']=dict(frequency=n,completed_quadrature_nodes=0,sums=sums)
            for lam in (-1,1):
                for j in range(Ntheta):
                    theta=2*math.pi*j/Ntheta;theta=(math.pi/2 if 4*j==Ntheta else 3*math.pi/2 if 4*j==3*Ntheta else theta)
                    row=expansion(n,lam,theta,M,J)
                    scale=math.pi*T*(1+lam/4)/Ntheta
                    for eps in p['epsilon']:
                        sd=sums[str(eps)]
                        paired=math.fsum(2*eps**k*moment[k]*row['coeff'][k] for k in range(2,len(row['coeff']),2))
                        closed=math.fsum(2*eps**k/(k+2)*row['coeff'][k] for k in range(2,len(row['coeff']),2))
                        for tau in gg:
                            b=eps*tau
                            plus=math.fsum(c*b**k for k,c in enumerate(row['coeff']))
                            minus=math.fsum(c*(-b)**k for k,c in enumerate(row['coeff']))
                            even=math.fsum(2*b**k*row['coeff'][k] for k in range(0,len(row['coeff']),2))
                            sol.require(abs(plus+minus-even)<=p['paired_polynomial_identity_max']*(1+abs(plus)+abs(minus)),'paired-sign polynomial identity',defect=abs(plus+minus-even))
                        src=math.fsum(2*eps**k*moment[k]*row['source_coeff_error'][k] for k in range(2,len(row['coeff']),2))
                        sd['correction']-=scale*paired
                        sd['source_force_error']+=abs(scale)*(T+rt)/T*src
                        sd['residual_force_error']+=abs(scale)*(T+rt)/T*row['remainder']
                        sd['primitive_error']+=abs(scale)*(T+rt)/T*row['primitive']
                        sd['path_polynomial_error']+=abs(scale)*abs(paired-closed)
                    report['partial_sums']['completed_quadrature_nodes']+=1;save()
                    if n==0:
                        subdivisions=p['angular_validation_subcells'];width=2*math.pi/Ntheta/subdivisions
                        for kk in range(subdivisions):
                            centre=theta-math.pi/Ntheta+(kk+.5)*width
                            expansion(n,lam,centre,M,J,angle_radius=width/2,validation=True)
            for eps in p['epsilon']:
                data=dict(sums[str(eps)])
                data['source_prefactor_error']=rt/(T-rt)*(abs(data['correction'])+data['source_force_error']+data['residual_force_error'])
                data['source_error']=data['source_force_error']+data['source_prefactor_error']
                data['primitive_error']+=floor*(1+abs(data['correction']))
                data['local_conditional_error']=math.fsum(data[k] for k in ('primitive_error','residual_force_error','path_polynomial_error'))
                data.update(n=n,epsilon=eps,configuration=[M,Ntheta,J,Nt],full_normalization=True,quadrature_angular_error_qualified=False,whole_frequency_sum=False)
                report['rows'].append(data)
            report['completed_frequencies'].append(dict(n=n,completed_quadrature_nodes=2*Ntheta,completed_expansion_nodes=nodes-before_nodes,continuous_validation_nodes=validationnodes-before_validation,calculation_cpu_seconds=time.process_time()-freq_cpu,wall_seconds=time.monotonic()-freq_start))
            report['partial_sums']={};save()
        rate0,rate32=report['completed_frequencies']
        report['empirical_finest_only_extrapolation']=dict(cpu_seconds=32*rate0['calculation_cpu_seconds']+480*rate32['calculation_cpu_seconds'],wall_seconds=32*rate0['wall_seconds']+480*rate32['wall_seconds'],certified=False,target_forecast_certified=False,assumption='All32 low-frequency finest slabs have n0 cost; all480 analytic finest slabs have n32 cost. Neither n dependence nor worst-case bound established. Excludes other hierarchy/mixed configurations and root launch/controller/intake costs.')
        report.update(outcome='complete_two_slab_supporting_pilot',pilot_complete=True,qualification_pass=False,whole_qualification_pass=False,target_release=False,failed_checks=[],next_step='Root reconciles genuine complete costs and independently reviews two-slab outcome. Empirical rate only; full qualification, original nonlinear C1/C2, whole forecast and parent coverage remain required before any target release.')
        guard()
    except Exception as error:
        report.update(outcome='partial_two_slab_obstacle',pilot_complete=False,qualification_pass=False,whole_qualification_pass=False,target_release=False,obstacle=dict(type=type(error).__name__,gate=getattr(error,'gate',str(error)),margins=getattr(error,'margins',{}),domain=context),failed_checks=['first_conditional_obstacle'],next_step='Root banks exact consumed partial slab coverage/costs and obtains independent materially different method/parent alternatives disposition; no replay or physical inference.')
    save();return report
