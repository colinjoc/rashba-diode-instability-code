"""Retrospective supporting reassembly, executable only under fresh native gates."""
import importlib.util,json,math,resource,sys,time,hashlib,copy
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
 sys.path.insert(0,str(package))
 p=json.loads((package/'QUALIFICATION_PROTOCOL.json').read_text());f=p['frozen'];v=[complex(a,b) for a,b in zip(f['v']['real'],f['v']['imag'])]
 eb=load(package/'error_bounds.py','response_bounds');ref=load(package/'source_reference.py','source_reference')
 out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.monotonic();rows=[];guard()
 binding=p['saved_qualification264'];path=package/binding['path']
 if hashlib.sha256(path.read_bytes()).hexdigest()!=binding['sha256']:raise ValueError('saved qualification hash')
 saved=json.loads(path.read_text());vn=sum(abs(x) for x in v)**2;reconstructed=saved['raw_force_reconstruction'];sb=p['saved_source256'];path=package/sb['path']
 if hashlib.sha256(path.read_bytes()).hexdigest()!=sb['sha256']:raise ValueError('saved source hash')
 sourcebudget=json.loads(path.read_text())
 if not sourcebudget['source_budget_pass'] or sourcebudget['integrated_source_budget']['Nw']!=p['levels'][1][-1]:raise ValueError('source bridge cutoff/pass')
 if sorted((r['Ntheta'],r['Nw']) for r in saved['rows'])!=sorted(map(tuple,p['schedule'])):raise ValueError('saved schedule')
 report=copy.deepcopy(saved);report.update(rows=rows,target_nodes=0,endpoint_confirmation_nodes=0,nonlinear_release=False,target_release=False,parent_answered=False,retrospective_reuse=True,new_response_evaluations=0,calculation_cpu_seconds=0,calculation_wall_seconds=0)
 report['schema']='project90105-retrospective-full-tail-qualification-v1'
 for old in saved['rows']:
  guard();row=copy.deepcopy(old);tail=eb.corrected_tail(f['x'],row['Nw'],f,p['source_box_radii'])
  row['tail']=tail;row['coefficient']=row['finite_coefficient']+tail['correction']
  # Original primitive already covered finite sums/old correction; add conservative new correction arithmetic separately.
  row['primitive_error']+=p['arithmetic_floor']*(1+abs(tail['correction']))
  row['numerical_error']=row['angular_error']+tail['remainder']+row['primitive_error'];rows.append(row)
 fine_tail=eb.corrected_tail(f['x'],p['levels'][1][-1],f,p['source_box_radii'])
 new_source_bound=sourcebudget['integrated_source_budget']['common_input_bound']+fine_tail['source_error']+p['arithmetic_floor']*(1+abs(fine_tail['correction']))
 report['source_bridge']={'original_common_bound':sourcebudget['integrated_source_budget']['common_input_bound'],'additional_full_tail_source_error':fine_tail['source_error'],'new_conservative_bound':new_source_bound,'domain_widened':False,'source_replayed':False,'reason':'Original mean-value enclosure retains finite4096 and old complete higher remainder variation. Add entire full C3 source variation, conservatively double-counting prior partial-C3 derivative allowance; old leading truncation variation also retained. No transfer to reference/current/K0.'}
 if True:
  by={(r['Ntheta'],r['Nw']):r for r in rows};A,B=p['levels'];fine=by[(A[-1],B[-1])];axes=[]
  for axis in (0,1):
   a=[by[(n,B[-1])] if axis==0 else by[(A[-1],n)] for n in p['levels'][axis]]
   axes.append(dict(axis=['theta','frequency'][axis],**hierarchy(*(r['coefficient'] for r in a),primitive=[r['primitive_error'] for r in a])))
  mixes=[];noises=[]
  for j in (0,1):
   keys=[(A[j],B[j]),(A[j],B[-1]),(A[-1],B[j]),(A[-1],B[-1])]
   mixes.append(by[keys[0]]['coefficient']-by[keys[1]]['coefficient']-by[keys[2]]['coefficient']+fine['coefficient']);noises.append(math.fsum(by[k]['primitive_error'] for k in keys))
  mixed=mixed_hierarchy(mixes,noises,[a['empirical_error'] for a in axes]);EI=fine['numerical_error']+math.fsum(a['empirical_error'] for a in axes)+mixed['empirical_error']
  reference=dict(saved['source_reference']);rawref=unpack(reference['raw_matrix']);correction,tailerror,td=ref.tail(f['x'],p['source_box_radii'],B[-1],f['K'])
  R=[[rawref[i][j]+correction[i][j] for j in range(2)] for i in range(2)]
  reference['matrix']=R;reference['raw_matrix']=rawref;reference['tail_entry_errors']=tailerror;reference['tail_details']=td
  reference['entry_errors']=[[reference['angular_entry_errors'][i][j]+tailerror[i][j]+reference['primitive_entry_errors'][i][j]+p['arithmetic_floor']*(1+abs(R[i][j])) for j in range(2)] for i in range(2)]
  raw=unpack(fine['raw_source_AP_matrix'])
  defects=[[abs(raw[i][j]-rawref[i][j]) for j in range(2)] for i in range(2)]
  discrepancy=max(ref.witness_error(defects,v),abs(contract(R,v)-fine['coefficient']));Eref=ref.witness_error(reference['entry_errors'],v)
  # No scalar current errors enter curvature. Independent operator reference
  # bounds plus measured mapping discrepancy + response EI form Ectrl.
  mapping_error=vn*fine['raw_source_to_force_mapping_discrepancy']
  Ectrl=max(p['historical_control_floor'],EI+Eref+discrepancy+mapping_error+p['arithmetic_floor']*(1+abs(contract(R,v))))
  report['source_reference']={k:(ref.pack(x) if k in ('matrix','raw_matrix') else x) for k,x in reference.items()};report['source_reference']['matrix_defects']=defects
  report['source_reference']['witness_error']=Eref;report['source_reference']['witness_discrepancy']=discrepancy
  report['normalization_controls']=saved['normalization_controls'];report['K0_controls']=saved['K0_controls']
  hermscale=p['hermiticity_relative_tolerance']*max(1.,max(abs(x) for row in raw for x in row));mapping=fine['raw_source_to_force_mapping_discrepancy']
  checks={'source256_coverage':new_source_bound<=f['Eep'],'signed_force_reconstruction':all(x['pass'] for x in reconstructed),'angular_hierarchy':axes[0]['pass'],'frequency_hierarchy':axes[1]['pass'],'mixed_hierarchy':mixed['pass'],'raw_Hermiticity':rawdefect(raw)<=hermscale,'raw_source_force_mapping':mapping<=p['normalization_tolerance'],'independent_reference_agreement':discrepancy<=fine['numerical_error']+Eref,'current_normalization':report['normalization_controls']['current_pass'],'mixed_normalization':report['normalization_controls']['mixed_pass'],'K0_amplitude_Ward':report['K0_controls']['pass'],'source_gap':report['K0_controls']['gap_selfconsistency_pass'],'control_adequacy':Ectrl<=min(f['EP']+f['Eep'],EI+f['Eep']),'primary_agreement':abs(fine['coefficient']-f['cP'])<=f['EP']+EI,'quarter_error':EI+f['Eep']<=.25*abs(fine['coefficient']),'negative_margin':fine['coefficient']+EI+f['Eep']<0,'original_primary_criteria':f['cP']+f['EP']+f['Eep']<0 and f['EP']+f['Eep']<=.25*abs(f['cP'])}
  report.update(axis_checks=axes,mixed_check=mixed,complete_physical_EI=EI,Ectrl=Ectrl,independent_coefficient=fine['coefficient'],checks=checks,failed_checks=[k for k,vv in checks.items() if not vv],qualification_pass=all(checks.values()),error_composition='EI=angular+tail+empirical primitive+axis+mixed; Ectrl=EI+independent entrywise reference witness error+measured paired-matrix discrepancy+measured force mapping contribution+mapping primitive. Eep retained once in acceptance, never in EI/Ectrl; scalar current errors excluded.')
  report['outcome']='complete_linear_qualification_pass_nonlinear_held' if report['qualification_pass'] else 'complete_linear_qualification_inadequate'
  report['obstacle']=None if report['qualification_pass'] else 'Named failed qualification conditions: '+', '.join(report['failed_checks'])
 report['calculation_cpu_seconds']=time.process_time()-cpu;report['calculation_wall_seconds']=time.monotonic()-wall;report['lifetime_rss_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 report['next_step']='Independent actual result review; linear pass only permits separately reviewed nonlinear C1/C2 design. Repeated unchanged failure requires independent parent alternatives; no physical answer.'
 raw=json.dumps(report,sort_keys=True,allow_nan=False).encode()+b'\n'
 if len(raw)>p['artifact_bytes']:raise RuntimeError('artifact ceiling')
 (out/'QUALIFICATION_RESULT.json').write_bytes(raw);return report
