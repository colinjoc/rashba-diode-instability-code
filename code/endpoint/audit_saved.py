"""Audit only saved observations. No new model evaluations or confirmation data."""
import json,hashlib,math,time
from pathlib import Path
import numpy as np

out=Path('science')
s=json.loads((out/'selected_checks.json').read_text());e=json.loads((out/'endpoint_certificate.json').read_text())
p=json.loads((out/'primary_evaluations.json').read_text());f=json.loads((out/'frozen_witness.json').read_text())
def lookup(x,K):
    rr=[r for r in p['records'] if np.array_equal(r['x'],x) and np.array_equal(r['K'],K)]
    assert rr
    return rr[-1]

slopes=[]
for axis in s['endpoint_propagation']['axes']:
    i=axis['i'];pr=axis['probes'];rho=e['rho'][i]
    b=[lookup(r['x'],s['K'])['table']['256,512']['b_matrix'] for r in pr]
    noise=(b[3]+b[0])/(2*rho)+(b[2]+b[1])/rho
    difference=abs(axis['full_centered']-axis['half_centered'])
    limit=.25*max(abs(axis['full_centered']),abs(axis['half_centered']))+noise
    slopes.append({'coordinate':i,'difference':difference,'rounding_only_noise':noise,'limit':limit,'pass':difference<=limit})

control=s['controls'];accepted_primary_error=s['EP']+s['Eep']
maximum_error=max(control['mixed_identity']['bound'],control['schur_identity']['bound'],control['amplitude_error'])
max_resolution_error=.25*abs(s['cP'])
rule=json.load(open('candidate/claim.json'))['protocol']['controls'][-1]
audit={'status':'INCONCLUSIVE: selected witness not accepted',
 'binding_sha256':'a21289d95d837d55cdef5e84a4eed1f752acf173248207cd802db8f05c251e6a',
 'audit_uses_saved_observations_only':True,'fixed_witness_sha256':hashlib.sha256((out/'frozen_witness.json').read_bytes()).hexdigest(),
 'preserved_history':{'selected_checks_json_original_primary_pass':s['primary_pass'],
     'explanation':'The saved preliminary flag omitted the final source-control error guard. It is retained as history and superseded by this audit. No nonlinear target perturbation was generated.'},
 'endpoint_propagation_smoothness':{'pass':all(a['pass'] for a in slopes),'axes':slopes,
     'correction':'Use only primitive-rounding bounds in the slope-contraction allowance. The original use of full EP was overly permissive; saved observations still pass the stricter check.'},
 'primary_sign':{'cP':s['cP'],'EP':s['EP'],'Eep':s['Eep'],'upper':s['primary_upper'],
                 'relative_error':s['relative_error'],'sign_and_25percent_guard_pass':True},
 'source_error_guard':{'frozen_rule':rule,'pass':bool(maximum_error<=accepted_primary_error),
     'accepted_primary_total_error':accepted_primary_error,
     'largest_error_allowed_by_primary_25percent_guard':max_resolution_error,
     'maximum_propagated_control_error':maximum_error,
     'mixed_control_error':control['mixed_identity']['bound'],
     'Schur_control_error':control['schur_identity']['bound'],
     'mixed_control_actual_discrepancy':control['mixed_identity']['difference'],
     'Schur_control_actual_discrepancy':control['schur_identity']['difference'],
     'interpretation':'Raw identity discrepancies are small. The conservative propagated numerical bounds are much larger than the admitted witness error. Counting only the discrepancy while discarding its prescribed bound would not supply the frozen numerical certificate.'},
 'accepted_instability':False,'independent_target_runs':0,
 'early_exit':'The implemented source-control certificate does not pass. Frozen early_exit classifies failed controls/incomplete certification as inconclusive; no alternative witness or expanded search is allowed.',
 'limits':'This audit does not prove that every equivalent implementation must produce the same loose bound, nor that the physical endpoint is stable or unstable. No independent curvature is available.',
 'hardware':{'probe':'science/cuda_probe_result.json','cuda_error':35,
              'message':'CUDA driver version is insufficient for CUDA runtime version',
              'nvml':'GPU access blocked by the operating system',
              'scope':'Additional execution limitation, not the cause of a scientific refutation. No claim about the external host driver is inferred.'},
 'followthrough':'Independent result review should assess this inconclusive certificate and the raw primary evidence. Any change to the source-control/error protocol requires prospective refinement and renewed review/admission. Preserve the fixed witness and observations; do not launch nonlinear confirmation or replace the witness under an assumption that controls passed.'}
assert not audit['source_error_guard']['pass']
(out/'certificate_audit.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n')
print(json.dumps(audit,indent=2))
