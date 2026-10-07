"""Fresh-admission guarded FINITE-BETA source/error qualification; target denied."""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
 os.environ[_name]='1'
import argparse,hashlib,importlib.util,json,resource,shutil,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);p.add_argument('--case',required=True);p.add_argument('--mode',choices=['qualification','target'],required=True);p.add_argument('--result',required=True);args=p.parse_args()
 if args.mode!='target':raise RuntimeError('This prospective package dispatches target mode only')
 work=Path(args.workspace).resolve();candidate=json.loads((work/'candidate/claim.json').read_text());binding=candidate['protocol']['execution_assets']['native_guard_adapter'];path=work/binding['path']
 if sha(path)!=binding['sha256']:raise RuntimeError('native guard binding mismatch')
 adapter=load(path,'native_linear_guard').Adapter(args);pp=adapter.asset('qualification_protocol');runner=adapter.asset('qualification_runner');definition=json.loads(pp.read_text());last=[-float('inf')]
 def full():
  adapter.guard();r=json.loads((adapter.root/'orchestrator/TARGET_EXECUTION_RELEASE.json').read_text())
  if r.get('scope')!='prospective_shifted_angular_fixed_witness_confirmation' or r.get('version')!=adapter.ctx['version'] or r.get('candidate_sha256')!=sha(adapter.candidate_path) or r.get('protocol_sha256')!=sha(pp) or r.get('worker_sha256')!=sha(Path(__file__)) or r.get('call_id')!=adapter.ctx['call_id']:raise RuntimeError('exact fresh finite-beta qualification release missing/changed')
  if shutil.disk_usage(work).free<definition['disk_reserve_bytes']:raise RuntimeError('disk reserve fails')
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>definition['local_memory_mb']*1024:raise RuntimeError('memory ceiling')
  last[0]=time.monotonic()
 def guard():
  if time.monotonic()>adapter.deadline-adapter.reserve:raise TimeoutError('finite allocation receipt reserve')
  if time.monotonic()-last[0]>=definition['full_guard_interval_seconds']:full()
 full()
 try:
  report=load(runner,'admitted_linear_qualification').run(pp.parent,adapter.out/'complete-finite-beta',guard);full()
  compact={k:report[k] for k in ['schema','outcome','qualification_pass','whole_qualification_pass','target_release','parent_answered','target_nodes','endpoint_confirmation_nodes','failed_checks','measurements']};compact['grid_count']=len(report['rows']);compact['details']='science/complete-finite-beta/QUALIFICATION_RESULT.json';compact['next_step']=report.get('next_step','Root owns genuine complete cost settlement and separate independent actual-result review; no automatic targetrelease.');compact['reused_293']=report.get('reused_293');compact['reused_299']=report.get('reused_299');compact['obstacle']=report.get('obstacle');adapter.finish(compact)
 except Exception as e:adapter.finish(None,e)
if __name__=='__main__':main()
