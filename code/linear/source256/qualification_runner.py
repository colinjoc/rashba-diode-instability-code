"""Isolated admitted common-source budget gate; no controls/response grids/targets."""
import importlib.util,json,resource,sys,time
from pathlib import Path
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def run(package,out,guard):
 p=json.loads((package/'QUALIFICATION_PROTOCOL.json').read_text());eb=load(package/'error_bounds.py','unchanged_source_tail');budget=load(package/'integrated_source_budget.py','admitted_integrated_source_budget')
 out.mkdir(exist_ok=False);cpu=time.process_time();wall=time.monotonic();f=p['frozen']
 report={'schema':'project90105-isolated-integrated-source-budget-v1','target_nodes':0,'endpoint_confirmation_nodes':0,'physical_response_evaluations':0,'rows':[],
         'qualification_pass':False,'source_budget_pass':False,'target_release':False,'parent_answered':False,'nonlinear_release':False,'obstacle':None,
         'outcome':'incomplete_isolated_source_budget','next_step':'Independent actual-result review; pass enables only a newly reviewed linear qualification protocol. Failure triggers independent material parent alternatives reassessment, not a bound/grid/cap retry.',
         'remaining_target_dependencies':p['remaining_target_dependencies']}
 try:
  guard();tail=eb.corrected_tail(f['x'],p['source_budget_frequency_cutoff'],f,p['source_box_radii']);guard()
  bound=budget.continuous_source_budget(f['x'],p['source_box_radii'],p['source_budget_frequency_cutoff'],p['source_budget_angular_cells'],f,tail,guard,p['arithmetic_floor']);guard()
  report['integrated_source_budget']=bound;report['frozen_shared_Eep']=f['Eep'];report['source_budget_pass']=bound['common_input_bound']<=f['Eep']
  report['qualification_pass']=report['source_budget_pass'] # Only this isolated gate; never linear or physical qualification.
  if not report['source_budget_pass']:report['obstacle']='Complete sufficient integrated source envelope exceeds unchanged Eep; no physical sign or unavoidable uncertainty established'
 except Exception as e:report['obstacle']=type(e).__name__+': '+str(e)
 report['outcome']='isolated_source_budget_pass_target_held' if report['source_budget_pass'] else 'isolated_source_budget_obstacle'
 report['calculation_cpu_seconds']=time.process_time()-cpu;report['calculation_wall_seconds']=time.monotonic()-wall;report['lifetime_rss_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 report['complete_physical_EI']=None;report['complete_Ectrl']=None;report['nonlinear_C1_C2']='NOT EXECUTED'
 raw=json.dumps(report,sort_keys=True,allow_nan=False).encode()+b'\n'
 if len(raw)>p['artifact_bytes']:raise RuntimeError('artifact ceiling')
 (out/'QUALIFICATION_RESULT.json').write_bytes(raw);return report
