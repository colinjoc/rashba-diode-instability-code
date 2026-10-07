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
    report=dict(schema='rashba-fixed-witness-origin-shift-confirmation-v1',scope=p['scope'],outcome='incomplete',qualification_pass=False,
        target_release=False,parent_answered=False,target_nodes=0,endpoint_confirmation_nodes=0,rows=[],failed_checks=[],
        data_role=p['data_role'],original_criterion=f,whole_qualification_pass=False,reused_293=[],new_control_primitive_overlap=True,primitive_empirical=True,directed_rounding=False,whole_memory_certified=False,
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
        publish_checkpoint(output/'QUALIFICATION_RESULT.json',raw,p['artifact_bytes'])
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
        savedpath=package/'SAVED293.json'
        sol.require(hashlib.sha256(savedpath.read_bytes()).hexdigest()==p['saved293_sha256'],'exact293 aggregate binding')
        saved=json.loads(savedpath.read_text())
        sol.require(saved['pilot_complete'] and saved['outcome']=='complete_two_slab_supporting_pilot' and not saved['failed_checks'],'genuine293 completed scoped controls')
        sol.require(saved['original_criterion']==f and len(saved['rows'])==6,'exact293 original source/witness/errors')
        sol.require(saved['accepted_linear_source_once']==p['accepted_linear_source_bound'] and saved['nonlinear_source_allowance']==p['nonlinear_source_allowance'],'exact293 source bridge reuse')
        reuse={}
        for sr in saved['rows']:
            sol.require(sr['configuration']==[16,64,8,8] and sr['n'] in [0,32] and sr['epsilon'] in p['epsilon'] and sr['full_normalization'],'exact293 observable domain')
            key=(sr['n'],sr['epsilon']);sol.require(key not in reuse,'unique293 observable');reuse[key]=sr
        sol.require(set(reuse)=={(n,e) for n in [0,32] for e in p['epsilon']},'complete293 exact slab reuse')
        originalpath=package/p['reuse293']['original_protocol']['path']
        sol.require(hashlib.sha256(originalpath.read_bytes()).hexdigest()==p['reuse293']['original_protocol']['sha256'],'exact293 protocol hash')
        original=json.loads(originalpath.read_text())
        for key in p['reuse293']['source_domain_keys']:
            sol.require(original[key]==p[key],'unchanged293 coverage domain',key=key)
        for name in ['nonlinear_response','correction_bounds']:
            sol.require(hashlib.sha256((package/(name+'.py')).read_bytes()).hexdigest()==p['prospective_'+name+'_sha256'],'exact prospective mathematics hash',name=name)
        sol.require(saved['completed_frequencies'][0]['n']==0 and saved['completed_frequencies'][0]['continuous_validation_nodes']==2048,'genuine293 continuous coverage count')
        report['reused_293']=dict(hash=p['saved293_sha256'],frequencies=[0,32],configuration=[16,64,8,8],rows=6,
            continuous_n0_validation_nodes=2048,original_costs_not_reset=True,new_confirmation=False,
            note='Nt8 totals authoritative. New Nt2/Nt4/subgrid/profile sufficient-statistic controls may share deterministic expansions; no pilot relaunch or discarded-coefficient fabrication.')
        # Completed299 rows are empirical historical evidence, not new confirmation.
        binding=p['reuse299']; saved299path=package/'SAVED299.json'
        sol.require(hashlib.sha256(saved299path.read_bytes()).hexdigest()==binding['report_sha256'],'exact failed299 saved checkpoint')
        priorprotocolpath=package/binding['original_protocol']['path']
        sol.require(hashlib.sha256(priorprotocolpath.read_bytes()).hexdigest()==binding['original_protocol']['sha256'],'exact frozen299 protocol')
        priorprotocol=json.loads(priorprotocolpath.read_text()); prior=json.loads(saved299path.read_text())
        for key in p['reuse293']['source_domain_keys']+['levels','schedule','mixed_pairs','finest','acceptance','nonlinear_source_allowance','accepted_linear_source_bound']:
            sol.require(priorprotocol[key]==p[key],'unchanged299 scientific definition',key=key)
        sol.require(prior['original_criterion']==f and not prior['failed_checks'],'299 completed controls have no saved scientific failures')
        sol.require(len(prior['rows'])==binding['completed_rows'],'complete retained299 row count')
        sol.require(prior['measurements']['continuous_validation_expansion_nodes']==binding['continuous_validation_nodes'],'retained299 continuous validation coverage')
        report['reused_299']=dict(continuous_validation_nodes=binding['continuous_validation_nodes'],scalar_rows_reused=0,full_circle=True,new_confirmation=False)
        report['reused_293']['scalar_rows_reused']=0
        report['reused_293']['rows']=0
        report['observation_identity']=p['observation_identity']
        baseline=json.loads((package/'SAVED308.json').read_text())
        sol.require(hashlib.sha256((package/'SAVED308.json').read_bytes()).hexdigest()==p['baseline308_sha256'],'qualification308 binding')
        sol.require(baseline['whole_qualification_pass'] and not baseline['failed_checks'],'complete308 baseline')
        # Genuine complete infinite-frequency/error bounds are checked first.
        tails={str(N):{str(e):bd.paired_tail(T,rt,D,rD,e,E,N,p['analytic_beta_radius'],floor) for e in p['epsilon']} for N in p['levels'][1]}
        report['complete_paired_infinite_tail_envelopes']=tails;report['analytic_tail_split']=p['analytic_tail_split'];save()
        mm={Nt:moments(Nt) for Nt in p['levels'][4]};gg={Nt:[float(t) for t in (np.polynomial.legendre.leggauss(Nt)[0]+1)/2] for Nt in p['levels'][4]};finest=p['finest']
        configs=[finest]+[c for c in p['schedule'] if c!=finest];base=[]
        for M,Nw,Ntheta,J,Nt in configs:
            key=(M,Ntheta,J)
            if key not in base:base.append(key)
        # Cache accumulated scalars only, never all coefficient arrays/history.
        for M,Ntheta,J in base:
            required=[c for c in configs if (c[0],c[2],c[3])==(M,Ntheta,J)]
            maxN=max(c[1] for c in required);Ntset=sorted(set(c[4] for c in required));sumdata={(eps,Nt):dict(correction=0.,source_force_error=0.,spatial_error=0.,primitive_error=0.,residual_force_error=0.,path_polynomial_error=0.,measured_tail_correction=0.) for eps in p['epsilon'] for Nt in Ntset}
            for n in range(maxN):
                activeNt=Ntset
                for lam in (-1,1):
                    for j in range(Ntheta):
                        theta=2*math.pi*j/Ntheta+p['angular_origin_radians']
                        if not activeNt:continue
                        report['current_evidence_role']='prospective_shifted_angular_confirmation'
                        row=expansion(n,lam,theta,M,J)
                        scale=math.pi*T*(1+lam/4)/Ntheta
                        for eps in p['epsilon']:
                            for Nt in activeNt:
                                sd=sumdata[(eps,Nt)];moment=mm[Nt]
                                # Pair cancellation is algebraic BEFORE absolute
                                # source/rounding/remainder aggregation; no odd term.
                                paired=math.fsum(2*eps**k*moment[k]*row['coeff'][k] for k in range(2,len(row['coeff']),2))
                                closed=math.fsum(2*eps**k/(k+2)*row['coeff'][k] for k in range(2,len(row['coeff']),2))
                                # Original paired-sign control, evaluated independently
                                # from the cancellation formula, on every Gauss node.
                                for tau in gg[Nt]:
                                    b=eps*float(tau)
                                    plus=math.fsum(c*b**k for k,c in enumerate(row['coeff']))
                                    minus=math.fsum(c*(-b)**k for k,c in enumerate(row['coeff']))
                                    even=math.fsum(2*b**k*row['coeff'][k] for k in range(0,len(row['coeff']),2))
                                    sol.require(abs(plus+minus-even)<=p['paired_polynomial_identity_max']*(1+abs(plus)+abs(minus)),'paired-sign polynomial identity',defect=abs(plus+minus-even))
                                sd['path_polynomial_error']+=abs(scale)*abs(paired-closed)
                                src=math.fsum(2*eps**k*moment[k]*row['source_coeff_error'][k] for k in range(2,len(row['coeff']),2))
                                sd['correction']-=scale*paired;sd['source_force_error']+=abs(scale)*(T+rt)/T*src
                                sd['residual_force_error']+=abs(scale)*(T+rt)/T*row['remainder']
                                sd['primitive_error']+=abs(scale)*(T+rt)/T*row['primitive']
                                if n>=p['analytic_tail_split']:sd['measured_tail_correction']-=scale*paired
                if (n+1) in [c[1] for c in required]:
                    for cfg in [c for c in required if c[1]==n+1]:
                        Nt=cfg[4]
                        for eps in p['epsilon']:
                            data=dict(sumdata[(eps,Nt)]);data['tail_error']=tails[str(n+1)][str(eps)]['bound']
                            data['source_prefactor_error']=rt/(T-rt)*(abs(data['correction'])+data['source_force_error']+data['residual_force_error']+data['tail_error'])
                            data['source_error']=data['source_force_error']+data['source_prefactor_error']
                            data['primitive_error']+=floor*(1+abs(data['correction']))
                            data['complete_numerical_error']=math.fsum(data[k] for k in ['spatial_error','primitive_error','residual_force_error','tail_error','path_polynomial_error'])
                            data.update(configuration=cfg,epsilon=eps,coefficient=anchor['independent_coefficient']+data['correction'],
                                complete_coefficient_error=anchor['complete_physical_EI']+data['complete_numerical_error'])
                            rows[(tuple(cfg),eps)]=data;report['rows'].append(data)
                    save()
        finest=p['finest'];all_checks=True;reductions=[];complete=[]
        for eps in p['epsilon']:
            axes=[]
            for axis,levels in enumerate(p['levels']):
                rr=[]
                for level in levels:cfg=list(finest);cfg[axis]=level;rr.append(rows[(tuple(cfg),eps)])
                x=hierarchy([r['coefficient'] for r in rr],[r['primitive_error'] for r in rr]);axes.append(x);all_checks &= x['pass_']
            mixed=[]
            for a,b in p['mixed_pairs']:
                cc=[];nn=[]
                for rung in (0,1):
                    ca=list(finest);cb=list(finest);cab=list(finest);ca[a]=cab[a]=p['levels'][a][rung];cb[b]=cab[b]=p['levels'][b][rung]
                    rr=[rows[(tuple(k),eps)] for k in [finest,ca,cb,cab]]
                    cc.append(rr[0]['coefficient']-rr[1]['coefficient']-rr[2]['coefficient']+rr[3]['coefficient']);nn.append(math.fsum(r['primitive_error'] for r in rr))
                ok=abs(cc[1])<=.75*abs(cc[0])+sum(nn) and abs(cc[1])<=.25*sum(x['error'] for x in axes)+nn[1]
                all_checks &= ok;mixed.append(dict(axes=[a,b],contrasts=cc,noise=nn,pass_=bool(ok),error=4*(abs(cc[1])+nn[1])))
            fine=rows[(tuple(finest),eps)];EI=fine['complete_coefficient_error']+sum(x['error'] for x in axes)+sum(x['error'] for x in mixed)
            source=fine['source_error'];budget=source<=p['nonlinear_source_allowance'];all_checks &= budget
            complete.append(dict(epsilon=eps,coefficient=fine['coefficient'],EI=EI,nonlinear_source_error=source,source_budget_pass=bool(budget)))
            reductions.append(dict(epsilon=eps,axes=axes,mixed=mixed))
        c=[x['coefficient'] for x in complete];b=[x['EI'] for x in complete];pr=[rows[(tuple(finest),e)]['primitive_error'] for e in p['epsilon']]
        d1,d2=c[1]-c[0],c[2]-c[1];noise=pr[0]+2*pr[1]+pr[2]
        epsok=abs(d2)<=.5*abs(d1)+noise and (d1*d2>=0 or min(abs(d1),abs(d2))<=noise)
        R01=(4*c[1]-c[0])/3;R12=(4*c[2]-c[1])/3
        er01=(4*b[1]+b[0])/3;er12=(4*b[2]+b[1])/3
        rich=4*(abs(R12-R01)+er01+er12)+er12
        # Original independent coefficient remains the separately qualified linear
        # anchor. Nonlinear C1/C2 are distinct complete negative-error gates.
        ci=anchor['independent_coefficient'];ei=anchor['complete_physical_EI'];ep=f['EP'];eep=f['Eep'];cp=f['cP'];ctrl=anchor['Ectrl']
        checks=dict(complete_schedule=len(report['rows'])==p['planned_rows'],epsilon_hierarchy=bool(epsok),
            original_linear_negative=ci+ei+eep<0 and cp+ep+eep<0,
            original_quarter_error=ei+eep<=.25*abs(ci) and ep+eep<=.25*abs(cp),
            original_primary_agreement=abs(cp-ci)<=ep+ei,
            original_control_adequacy=ctrl<=min(ep+eep,ei+eep),
            all_source_controls=bool(anchor['qualification_pass']),
            complete_nonlinear_C1_negative=c[1]+b[1]+eep<0,
            complete_nonlinear_C2_negative=c[2]+b[2]+eep<0,
            complete_richardson_negative=R12+rich+eep<0,
            source_budget=all(x['source_budget_pass'] for x in complete),
            complete_axis_mixed=bool(all_checks))
        origin_comparisons=[]
        for fresh,old in zip(complete,baseline['complete_results']):
            sol.require(fresh['epsilon']==old['epsilon'],'exact308 comparison epsilon')
            margin=fresh['EI']+old['EI']
            origin_comparisons.append(dict(epsilon=fresh['epsilon'],difference=abs(fresh['coefficient']-old['coefficient']),allowed=margin,pass_=abs(fresh['coefficient']-old['coefficient'])<=margin))
        checks['shifted_origin_agreement']=all(x['pass_'] for x in origin_comparisons)
        report['origin_comparisons']=origin_comparisons
        report['endpoint_confirmation_nodes']=nodes
        report['target_nodes']=nodes
        report.update(complete_results=complete,hierarchy=reductions,epsilon_reduction=dict(C1=c[1],C2=c[2],Richardson01=R01,Richardson12=R12,complete_Richardson_error=rich,pass_=bool(epsok)),
            EI=ei,Ectrl=ctrl,checks=checks,failed_checks=[k for k,v in checks.items() if not v],
            qualification_pass=all(checks.values()),whole_qualification_pass=all(checks.values()),outcome='supported_empirical_fixed_witness' if all(checks.values()) else 'inconclusive_target',
            next_step='Root owns genuine accounting, independent target result and parent coverage review; no automatic parent closure. Inadequacy returns exact evidence for materially different B/C strategy review.')
        report['complete_test_forecast']=dict(complete_schedule_measured=True,calculation_cpu_seconds=time.process_time()-cpu,
            calculation_wall_seconds=time.monotonic()-start,all_original_rows=len(report['rows'])==p['planned_rows'],
            pending_root_settlement=True,target_forecast_certified=False,target_release=False,
            scope='Actual complete supporting qualification only. Outer launch/controller/intake and independent prospective endpoint protocol forecast remain root-owned; no endpoint release from this field.')
        guard()
    except Exception as error:
        report.update(outcome='precise_conditional_inclusion_or_resource_obstacle',obstacle=dict(type=type(error).__name__,gate=getattr(error,'gate',str(error)),margins=getattr(error,'margins',{}),domain=context),
            failed_checks=['first_conditional_obstacle'],next_step='Root banks consumed finite-beta B failure, exact domain and costs; independent material B/C disposition. No cosmetic retry, target or physical inference.')
    save();return report
