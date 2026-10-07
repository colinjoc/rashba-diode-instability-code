"""Complete prospective fixed-witness linear qualification, target withheld."""
import importlib.util,json,math,resource,sys,time,hashlib
from pathlib import Path

def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def hierarchy(a,b,c,primitive):
 b0,b1,b2=primitive;d1=b-a;d2=c-b;noise=b0+2*b1+b2
 return {'d1':d1,'d2':d2,'noise':noise,'pass':abs(d2)<=.75*abs(d1)+noise and (d1*d2>=0 or min(abs(d1),abs(d2))<=noise),'empirical_error':4*(abs(d2)+b1+b2)+b2}

def mixed_hierarchy(contrasts,primitive_noises,axis_errors):
 lo,hi=contrasts;blo,bhi=primitive_noises
 return {'contrasts':contrasts,'primitive_noises':primitive_noises,'pass':abs(hi)<=.75*abs(lo)+blo+bhi and abs(hi)<=.25*math.fsum(axis_errors)+bhi,'empirical_error':4*(abs(hi)+bhi)}

def unpack(M):return [[complex(*x) for x in row] for row in M]
def rawdefect(M):return max(abs(M[i][j]-M[j][i].conjugate()) for i in range(2) for j in range(2))
def contract(M,v):return math.fsum((v[i].conjugate()*M[i][j]*v[j]).real for i in range(2) for j in range(2))

def run(package,out,guard):
 p=json.loads((package/'QUALIFICATION_PROTOCOL.json').read_text());f=p['frozen'];T,D,q,h=f['x'];v=[complex(a,b) for a,b in zip(f['v']['real'],f['v']['imag'])]
 lin=load(package/'linear_response.py','response_solver');eb=load(package/'error_bounds.py','response_bounds');ref=load(package/'source_reference.py','source_reference');cur=load(package/'current_controls.py','current_controls')
 out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.monotonic();rows=[]
 report={'schema':'project90105-allK-fixed-witness-linear-qualification-v1','rows':rows,'qualification_pass':False,'target_release':False,'nonlinear_release':False,'parent_answered':False,'target_nodes':0,'endpoint_confirmation_nodes':0,'obstacle':None,'nonlinear_C1_C2':'NOT EXECUTED: mandatory fresh independent protocol/admission/root release','remaining_target_dependencies':p['remaining_target_dependencies']}
 try:
  guard();binding=p['saved_source256'];savedpath=package/binding['path']
  if hashlib.sha256(savedpath.read_bytes()).hexdigest()!=binding['sha256']:raise ValueError('saved256 hash mismatch')
  saved=json.loads(savedpath.read_text());budget=saved['integrated_source_budget'];report['saved_source256']={'sha256':binding['sha256'],'common_input_bound':budget['common_input_bound'],'replayed':False,'coverage':'Original fixed witness only, unchanged full source box; cutoff4096/source C3 and remainder semantics preserved'}
  if not saved['source_budget_pass'] or budget['Nw']!=p['levels'][1][-1] or budget['common_input_bound']>f['Eep']:raise ValueError('saved source256 coverage mismatch')
  source=lin.SourceParameters(T,h,q,D,tuple(f['K']));rt,rD,rq,rh=p['source_box_radii'];tm=T-rt;dm=(abs(q)+rq)/2+abs(h)+rh;km=math.hypot(*f['K']);vn=sum(abs(x) for x in v)**2
  reconstructed=[]
  for helicity in (-1,1):
   for c in (1.,0.,-1.):
    z=complex(math.pi*T,(q/2-helicity*h)*c);projection=f['K'][0]*c
    B,H=lin.force_matrix(z,D,projection);plus=lin.harmonic_response(z,D,projection,1);minus=lin.harmonic_response(z,D,projection,-1)
    ep=v[0]+1j*v[1];em=v[0].conjugate()+1j*v[1].conjugate()
    force=(ep.conjugate()*(plus['f_a']*v[0]+plus['f_p']*v[1])+em.conjugate()*(minus['f_a']*v[0].conjugate()+minus['f_p']*v[1].conjugate())).real
    expected=lin.contract(H,v);primitive=p['arithmetic_floor']*(1+abs(force)+abs(expected))
    reconstructed.append({'helicity':helicity,'cos_theta':c,'force':force,'H_contraction':expected,'discrepancy':abs(force-expected),'primitive_error':primitive,'raw_B':ref.pack(B),'pass':abs(force-expected)<=p['normalization_tolerance']+primitive,'scope':'Actual signed Fourier force reconstruction from both sectors; raw B is not Hermitian kernel'})
  report['raw_force_reconstruction']=reconstructed
  for Nt,Nw in p['schedule']:
   guard();r=lin.finite_grid_diagnostic(source,angular_nodes=Nt,frequencies=Nw,guard=guard,witness=v,max_evaluations=p['max_evaluations_per_grid'])
   angular=math.pi*(T+rt)*2*math.fsum(eb.angular_strip_bound((2*n+1)*math.pi*tm,dm,D-rD,D+rD,km,Nt,vn)['bound'] for n in range(Nw));tail=eb.corrected_tail(f['x'],Nw,f,p['source_box_radii'])
   # Primitive floor uses own accumulated absolute-entry majorant; no sourcebox replay.
   absolute=math.pi*(T+rt)*2*vn*math.fsum(2/((2*n+1)*math.pi*tm)+eb.derivative_bounds((2*n+1)*math.pi*tm,dm,D+rD,km,D-rD)['H'] for n in range(Nw))
   primitive=p['arithmetic_floor']*(1+absolute+abs(2*math.log(tm))*vn+abs(tail['correction']))
   c=r['finite_grid_witness_coefficient']+tail['correction'];EI=angular+tail['remainder']+primitive
   rows.append({'Ntheta':Nt,'Nw':Nw,'coefficient':c,'finite_coefficient':r['finite_grid_witness_coefficient'],'finite_matrix':r['matrix'],'raw_source_AP_matrix':r['raw_source_AP_matrix'],'raw_Hermiticity_defect':r['raw_source_Hermiticity_defect'],'raw_source_to_force_mapping_discrepancy':r['raw_source_to_force_mapping_discrepancy'],'evaluations':r['grid']['evaluations'],'angular_error':angular,'tail':tail,'primitive_error':primitive,'numerical_error':EI})
  by={(r['Ntheta'],r['Nw']):r for r in rows};A,B=p['levels'];fine=by[(A[-1],B[-1])];axes=[]
  for axis in (0,1):
   a=[by[(n,B[-1])] if axis==0 else by[(A[-1],n)] for n in p['levels'][axis]]
   axes.append(dict(axis=['theta','frequency'][axis],**hierarchy(*(r['coefficient'] for r in a),primitive=[r['primitive_error'] for r in a])))
  mixes=[];noises=[]
  for j in (0,1):
   keys=[(A[j],B[j]),(A[j],B[-1]),(A[-1],B[j]),(A[-1],B[-1])]
   mixes.append(by[keys[0]]['coefficient']-by[keys[1]]['coefficient']-by[keys[2]]['coefficient']+fine['coefficient']);noises.append(math.fsum(by[k]['primitive_error'] for k in keys))
  mixed=mixed_hierarchy(mixes,noises,[a['empirical_error'] for a in axes]);EI=fine['numerical_error']+math.fsum(a['empirical_error'] for a in axes)+mixed['empirical_error']
  reference=ref.evaluate(f['x'],f['K'],A[-1],B[-1],p['source_box_radii'],guard,p['arithmetic_floor']);R=reference['matrix'];raw=unpack(fine['raw_source_AP_matrix']);rawref=reference['raw_matrix']
  defects=[[abs(raw[i][j]-rawref[i][j]) for j in range(2)] for i in range(2)]
  discrepancy=ref.witness_error(defects,v);Eref=ref.witness_error(reference['entry_errors'],v)
  # No scalar current errors enter curvature. Independent operator reference
  # bounds plus measured mapping discrepancy + response EI form Ectrl.
  mapping_error=vn*fine['raw_source_to_force_mapping_discrepancy']
  Ectrl=max(p['historical_control_floor'],EI+Eref+discrepancy+mapping_error+p['arithmetic_floor']*(1+abs(contract(R,v))))
  report['source_reference']={k:(ref.pack(x) if k in ('matrix','raw_matrix') else x) for k,x in reference.items()};report['source_reference']['matrix_defects']=defects
  report['source_reference']['witness_error']=Eref;report['source_reference']['witness_discrepancy']=discrepancy
  report['normalization_controls']=cur.evaluate(p,ref,guard)
  # Actual nonzero-source K0 amplitude/gap/Ward independent code comparison.
  zero=lin.finite_grid_diagnostic(lin.SourceParameters(T,h,q,D,(0.,0.)),angular_nodes=A[-1],frequencies=B[-1],guard=guard,max_evaluations=p['max_evaluations_per_grid']);zr=ref.evaluate(f['x'],[0.,0.],A[-1],B[-1],p['source_box_radii'],guard,p['arithmetic_floor'])
  zd=[[abs(unpack(zero['raw_source_AP_matrix'])[i][j]-zr['raw_matrix'][i][j]) for j in range(2)] for i in range(2)]
  # K0 derivative identities are direct source-functional G definitions,
  # independently bound and compared to both Riccati sectors numerically.
  controlscale=p['source_control_relative_tolerance']*max(1.,abs(zr['matrix'][0][0]),abs(zr['matrix'][1][1]));zerr=max(max(row) for row in zr['entry_errors'])
  report['K0_controls']={'amplitude_FDD':zr['matrix'][0][0].real,'phase_FD_over_D':zr['matrix'][1][1].real,'gap_residual':D*zr['matrix'][1][1].real/2,'gap_residual_error':D*zr['entry_errors'][1][1]/2,'entry_errors':zr['entry_errors'],'finite_regulator_defects':zd,'raw_paired_matrix':zero['raw_source_AP_matrix'],'pass':max(max(row) for row in zd)<=p['normalization_tolerance']+2*zerr and zerr<=controlscale,'identity':'M_AA=F_DD; M_PP=F_D/D; gap residual preserved, no forced Ward zero','gap_selfconsistency_pass':abs(D*zr['matrix'][1][1].real/2)<=controlscale+D*zr['entry_errors'][1][1]/2}
  hermscale=p['hermiticity_relative_tolerance']*max(1.,max(abs(x) for row in raw for x in row));mapping=fine['raw_source_to_force_mapping_discrepancy']
  checks={'source256_coverage':True,'signed_force_reconstruction':all(x['pass'] for x in reconstructed),'angular_hierarchy':axes[0]['pass'],'frequency_hierarchy':axes[1]['pass'],'mixed_hierarchy':mixed['pass'],'raw_Hermiticity':rawdefect(raw)<=hermscale,'raw_source_force_mapping':mapping<=p['normalization_tolerance'],'independent_reference_agreement':discrepancy<=fine['numerical_error']+Eref,'current_normalization':report['normalization_controls']['current_pass'],'mixed_normalization':report['normalization_controls']['mixed_pass'],'K0_amplitude_Ward':report['K0_controls']['pass'],'source_gap':report['K0_controls']['gap_selfconsistency_pass'],'control_adequacy':Ectrl<=min(f['EP']+f['Eep'],EI+f['Eep']),'primary_agreement':abs(fine['coefficient']-f['cP'])<=f['EP']+EI,'quarter_error':EI+f['Eep']<=.25*abs(fine['coefficient']),'negative_margin':fine['coefficient']+EI+f['Eep']<0,'original_primary_criteria':f['cP']+f['EP']+f['Eep']<0 and f['EP']+f['Eep']<=.25*abs(f['cP'])}
  report.update(axis_checks=axes,mixed_check=mixed,complete_physical_EI=EI,Ectrl=Ectrl,independent_coefficient=fine['coefficient'],checks=checks,failed_checks=[k for k,vv in checks.items() if not vv],qualification_pass=all(checks.values()),error_composition='EI=angular+tail+empirical primitive+axis+mixed; Ectrl=EI+independent entrywise reference witness error+measured paired-matrix discrepancy+measured force mapping contribution+mapping primitive. Eep retained once in acceptance, never in EI/Ectrl; scalar current errors excluded.')
  report['outcome']='complete_linear_qualification_pass_nonlinear_held' if report['qualification_pass'] else 'complete_linear_qualification_inadequate'
  if not report['qualification_pass']:report['obstacle']='Named failed qualification conditions: '+', '.join(report['failed_checks'])
 except Exception as e:
  report['obstacle']=type(e).__name__+': '+str(e);report['outcome']='linear_qualification_execution_obstacle';report['complete_error_status']='Not calculated because exact named execution obstacle prevented complete observable; no sign inference'
 report['calculation_cpu_seconds']=time.process_time()-cpu;report['calculation_wall_seconds']=time.monotonic()-wall;report['lifetime_rss_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 report['next_step']='Independent result review; pass enables only separately reviewed nonlinear C1/C2 protocol; inadequacy needs one materially different evidence-bound remedy or parent alternatives disposition. No physical parent answer.'
 raw=json.dumps(report,sort_keys=True,allow_nan=False).encode()+b'\n'
 if len(raw)>p['artifact_bytes']:raise RuntimeError('artifact ceiling')
 (out/'QUALIFICATION_RESULT.json').write_bytes(raw);return report
