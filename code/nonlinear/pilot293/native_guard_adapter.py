#!/usr/bin/env python3
"""Inert-on-import deterministic adapter for genuine external engine science.

Only the root controller dispatches after exact fresh admission. Qualification
never dispatches targets. A target chunk requires a same-version real reviewed
release. No case, ledger, admission or release writes occur here.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time


def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def read(p): return json.loads(Path(p).read_text())

_STATE_CACHE = None

def _state_fingerprint(stat):
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)

def cached_state(path):
    """Trust local mutation metadata; stable before/after or fail closed.

    Only JSON parsing is cached. Call/version/admission, exact candidate/asset
    hashes, authority/STOP, deadlines and CPU checks are still freshly run.
    Root owns state writes; cached dictionaries are read-only at these sites.
    """
    global _STATE_CACHE
    path = Path(path)
    before = _state_fingerprint(path.stat())
    if _STATE_CACHE is not None and _STATE_CACHE[0] == str(path) and _STATE_CACHE[1] == before:
        if _state_fingerprint(path.stat()) != before:
            _STATE_CACHE = None
            raise RuntimeError('state changed during cached read')
        return _STATE_CACHE[2]
    _STATE_CACHE = None
    with path.open('rb') as stream:
        if _state_fingerprint(os.fstat(stream.fileno())) != before:
            raise RuntimeError('state changed before parse')
        raw = stream.read()
        if _state_fingerprint(os.fstat(stream.fileno())) != before:
            raise RuntimeError('state changed during read')
    value = json.loads(raw)
    if _state_fingerprint(path.stat()) != before:
        raise RuntimeError('state changed during parse')
    _STATE_CACHE = (str(path), before, value)
    return value




def cpu():
    a=resource.getrusage(resource.RUSAGE_SELF); b=resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime+a.ru_stime+b.ru_utime+b.ru_stime


def outer_cpu_sample(ctx, active):
    """Actual bound controller SELF and awaited-other CHILD ticks, while worker lives."""
    pid=int(ctx['controller_pid'])
    if pid!=int(active['pid']) or str(ctx['controller_process_start'])!=str(active['controller_process_start']):
        raise RuntimeError('outer controller identity differs from actual allocation')
    row=(Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
    if row[19]!=str(ctx['controller_process_start']):
        raise RuntimeError('outer controller PID/start identity changed')
    ticks=int(os.sysconf('SC_CLK_TCK'))
    if ticks<=0:raise RuntimeError('invalid CPU tick frequency')
    values=[int(row[i]) for i in (11,12,13,14)]
    if min(values)<0:raise RuntimeError('invalid outer CPU counters')
    return dict(self_seconds=(values[0]+values[1])/ticks,
                awaited_other_children_seconds=(values[2]+values[3])/ticks,
                quantization_allowance_seconds=2./ticks)


def outer_cpu_delta(now, before):
    parts=[]
    for key in ('self_seconds','awaited_other_children_seconds'):
        delta=finite(now[key])-finite(before[key])
        allowance=finite(now['quantization_allowance_seconds'])
        if delta < -allowance:raise RuntimeError('outer CPU counter decreased')
        parts.append(max(0.,delta)+allowance)
    return sum(parts)


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def new_json(path,data,cap=65536):
    raw=(json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    if len(raw)>cap: raise ValueError('64KiB record ceiling exceeded')
    with open(path,'xb') as f:f.write(raw)
    return dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def finite(x):
    x=float(x)
    if not math.isfinite(x) or x<0: raise ValueError('invalid nonnegative certified bound')
    return x


class Adapter:
    def __init__(self,args):
        self.args=args; self.work=Path(args.workspace).resolve(); self.root=Path(args.case).resolve()
        self.out=self.work/'science';self.out.mkdir(exist_ok=True)
        self.start=time.monotonic();self.cpu_start=cpu()
        self.candidate_path=self.work/'candidate/claim.json';self.candidate=read(self.candidate_path)
        self.protocol=self.candidate['protocol'];self.ctx=read(self.work/'execution_context.json')
        self.assets=self.protocol['execution_assets']; self.asset_paths={}
        self.inputs=cached_state(self.root/'state.json')['calls'][int(self.ctx['call_id'])-1]['inputs']
        # The parent runtime's deadline can precede the engine deadline.
        ledger=read(self.root/'orchestrator/computational-budget.json');active=ledger['active_allocation']
        if active is None or active['version']!=self.ctx['version'] or active['next_call_id']!=self.ctx['call_id']:
            raise RuntimeError('actual root-owned active allocation missing or changed')
        self.deadline=min(float(self.ctx['deadline_monotonic']),float(active['whole_allocation_deadline_monotonic']))
        self.reserve=90 if args.mode=='target' else 45
        self.baseline_ledger=ledger; self.active=active
        self.outer_allocation_baseline=dict(
            self_seconds=finite(active['controller_self_cpu_at_allocation_start']),
            awaited_other_children_seconds=finite(active['controller_children_cpu_at_allocation_start']))
        outer_cpu_sample(self.ctx,self.active) # missing/changed identity/counters fails before science
        self.new=[];self.new_summaries=[];self.terminal=False
        repo=Path(self.root).parents[2];sys.path.insert(0,str(repo))
        from engine.case_authority import require_authority
        from orchestration.physics100.novelty import require_admission
        self.require_authority=require_authority;self.require_admission=require_admission
        self.guard()

    def asset(self,name):
        b=self.assets[name];rel=Path(b['path'])
        if rel.is_absolute():
            if name!='nvcc': raise RuntimeError('only exact hashed nvcc tooling may be absolute')
            path=rel
        else:
            if '..' in rel.parts:raise RuntimeError('asset escapes immutable workspace')
            path=self.work/rel
            if self.inputs.get(str(rel))!=b['sha256']:raise RuntimeError('asset absent from actual call input binding: '+name)
        if not path.is_file() or sha(path)!=b['sha256']:raise RuntimeError('asset missing/changed: '+name)
        self.asset_paths[name]=path;return path

    def current_ledger(self):
        l=read(self.root/'orchestrator/computational-budget.json')
        a=l.get('active_allocation')
        if not a or a['version']!=self.ctx['version'] or a['next_call_id']!=self.ctx['call_id'] or l.get('budget_violation'):
            raise RuntimeError('controller allocation no longer valid')
        return l

    def guard(self,phase=None,record=None):
        self.require_authority(self.root,self.root.name)
        s=cached_state(self.root/'state.json')
        if s['active_call']!=self.ctx['call_id'] or s['current_version']!=self.ctx['version']:
            raise RuntimeError('engine call/version changed')
        v=s['versions'][int(self.ctx['version'])-1]
        if self.ctx['backend']!='external' or self.ctx['stage']!='science' or s['admission'] is None:
            raise RuntimeError('fresh admitted genuine external science required')
        if sha(self.root/v['path']/'claim.json')!=sha(self.candidate_path):raise RuntimeError('candidate changed')
        self.require_admission(v['project'],self.root/v['path']/'admission')
        l=self.current_ledger()
        if bool(l['active_allocation']['qualification_only'])!=(self.args.mode=='qualification'):
            raise RuntimeError('mode differs from actual controller allocation')
        if time.monotonic()>self.deadline-self.reserve:raise TimeoutError('receipt/intake reserve reached')
        if cpu()-self.cpu_start+outer_cpu_delta(outer_cpu_sample(self.ctx,self.active),self.outer_allocation_baseline)>=float(l['active_allocation']['per_invocation_cpu_cap_seconds']):
            raise RuntimeError('local CPU allowance exhausted')
        for name in self.assets:self.asset(name)
        if self.args.mode=='target':self.check_release()

    def check_release(self):
        r=read(self.root/'orchestrator/TARGET_EXECUTION_RELEASE.json')
        if r['version']!=self.ctx['version'] or r['candidate_sha256']!=sha(self.candidate_path):
            raise RuntimeError('stale target release')
        s=cached_state(self.root/'state.json');review=s['calls'][int(r['independent_review_call'])-1]
        path=self.root/'calls'/f"{review['id']:06d}"/'receipt.json'
        if review['stage']!='result_review' or review['version']!=self.ctx['version']:
            raise RuntimeError('target release lacks same-version actual result review')
        if sha(path)!=r['independent_review_receipt_sha256'] or read(path)['result']['decision']!='continue':
            raise RuntimeError('target release continuation binding invalid')
        return r

    def allowances(self):
        # Constants are prospective reviewed bounds, NEVER measured claims.
        a=self.protocol['resource_projection'].get('bounded_overhead_allowances')
        expected=dict(per_node_host_cpu_seconds=1.2,remaining_main_cpu=100.,remaining_diagnostic_cpu=50.,
                      remaining_other_cpu=250.,remaining_gpu_seconds=100.,remaining_wall_seconds=1500.,
                      checkpoint_wall_seconds=660.)
        if not isinstance(a,dict) or any(float(a.get(k,-1))!=v for k,v in expected.items()):
            raise RuntimeError('exact reviewed nonzero overhead allowances missing')
        if 896*a['per_node_host_cpu_seconds']+a['remaining_main_cpu']+a['remaining_diagnostic_cpu']+a['remaining_other_cpu']>1500:
            raise RuntimeError('1500CPU overhead envelope violated')
        return a

    def qualification(self):
        self.guard();a=self.allowances();cert=read(self.asset('controller_certificate'))
        if shutil.disk_usage(self.work).free<256*1024*1024:
            raise RuntimeError('disk free below256MiB; no cleanup or science rerun authorized')
        ledger=self.current_ledger();remaining_setup=finite(cert['remaining_original_setup_cpu_seconds'])
        # Snapshot grants no replenishment: named later qualification/preparation debits required.
        debit=finite(self.protocol['resource_projection']['setup_debits_after_certificate_cpu_seconds'])
        remaining_setup-=debit
        if remaining_setup<=0:raise RuntimeError('remaining original setup reserve exhausted')
        role_names=('direct_cuda_source','direct_root_fixture_report','direct_source_audit','v6_cuda_source','nvcc','reference_generic_raw','reference_generic_audit',
                    'reference_stress_binary','driver','reducer','sums_library','controller_certificate')
        roles={n:dict(path=str(self.asset(n)),sha256=self.assets[n]['sha256']) for n in role_names}
        roles['candidate']=dict(path=str(self.candidate_path),sha256=sha(self.candidate_path))
        roles['execution_context']=dict(path=str(self.work/'execution_context.json'),sha256=sha(self.work/'execution_context.json'))
        # An admitted exact candidate bounds allowances; bind it rather than invent an unreviewed certificate.
        overhead=dict(controller_certified=True,certificate_binding=roles['candidate'],
                      remaining_main_cpu=a['remaining_main_cpu'],remaining_diagnostic_cpu=a['remaining_diagnostic_cpu'],
                      remaining_other_cpu=a['remaining_other_cpu'],remaining_gpu_seconds=a['remaining_gpu_seconds'],
                      remaining_wall_seconds=a['remaining_wall_seconds'],checkpoint_wall_seconds=a['checkpoint_wall_seconds'],
                      current_generated_bytes=finite(cert['existing_generated_bytes']),
                      remaining_generated_bytes_bound=finite(cert['remaining_generated_bytes_allowance']))
        envelope=ledger['prospective_device_wall']
        manifest=dict(repository_root=str(self.root.parents[2]),roles=roles,scope='qualification_only',
                      certified_remaining_setup_cpu_seconds=remaining_setup,
                      charged_cpu_since_certificate=max(0.,ledger['charged_cpu_seconds']-cert['charged_cpu_seconds'])+(cpu()-self.cpu_start),
                      new_release_science_wall_charged_seconds=envelope['charged_science_wall_seconds'],
                      per_node_host_cpu_allowance_seconds=a['per_node_host_cpu_seconds'],overhead_bounds=overhead)
        folder=self.out/'qualification';folder.mkdir(exist_ok=True)
        mf=self.out/'qualification_worker_manifest.json';new_json(mf,manifest)
        worker=load_module(self.asset('qualification_worker'),'actual_qualification_worker')
        argv=sys.argv
        try:
            sys.argv=[str(self.asset('qualification_worker')),'--case-root',str(self.root),'--workspace',str(folder),'--manifest',str(mf)]
            code=worker.main()
        finally:sys.argv=argv
        report=read(folder/'RESULT.json')
        if report.get('target_nodes')!=0:raise RuntimeError('qualification generated forbidden targets')
        # A reported FAIL is a valid scientific result, not a transport failure.
        return dict(mode='qualification',outcome='supporting_qualification_pass' if report['qualification_pass'] else 'inconclusive',
                    qualification_pass=report['qualification_pass'],target_nodes=0,
                    report=report,worker_exit_code=code,report_sha256=sha(folder/'RESULT.json'),
                    next_step='Actual independent result review and exact root target release required; qualification never admits target science.')

    def prior_file(self,call_id,digest):
        base=self.work/'prior_artifacts'/f'{int(call_id):06d}'
        matches=[]
        for path in base.rglob('*'):
            if path.is_file() and sha(path)==digest:
                rel=str(path.relative_to(self.work))
                if self.inputs.get(rel)!=digest:raise RuntimeError('unbound prior artifact')
                matches.append(path)
        if len(matches)!=1:raise RuntimeError('missing or ambiguous actual qualification artifact '+digest)
        return matches[0]

    def validated_prior_nodes(self,expected,source_hash,binary_hash):
        found={};qualified_hashes={source_hash,binary_hash,self.assets['driver']['sha256'],
                                  self.assets['reducer']['sha256'],self.assets['sums_library']['sha256']}
        for path in sorted((self.work/'prior_artifacts').rglob('*.json')):
            # Reserved target-summary paths are evidence even if their JSON or
            # schema is corrupt. Never silently drop them and replay a node.
            reserved=(path.parent.name=='nodes' and path.parent.parent.name=='science'
                      and not path.name.endswith('.execution.json'))
            try:s=read(path)
            except (ValueError,OSError) as error:
                if reserved:raise RuntimeError('malformed recognized target summary: '+str(path)) from error
                continue
            if not isinstance(s,dict):
                if reserved:raise RuntimeError('recognized target summary is not an object: '+str(path))
                continue
            if s.get('schema')!='proposed-compact-rashba-node-v1':
                if reserved:raise RuntimeError('recognized target summary schema missing/changed: '+str(path))
                continue
            ident=s.get('identity',{})
            if not isinstance(ident,dict):raise RuntimeError('recognized target identity is not an object')
            key=ident.get('id')
            if ident.get('family')=='qualification':
                if reserved:raise RuntimeError('qualification cannot occupy reserved target-summary path')
                continue
            if key not in expected or ident!=expected[key]:raise RuntimeError('unknown or altered target node')
            rel=str(path.relative_to(self.work));digest=sha(path)
            if self.inputs.get(rel)!=digest:raise RuntimeError('target summary input hash missing')
            if key in found:raise RuntimeError('duplicate immutable target node identity')
            if s.get('valid') is not True or s.get('payload_bytes')!=7864320 or s.get('record_count')!=65536:
                raise RuntimeError('invalid saved compact target record')
            if type(s.get('maxima',{}).get('passes')) is not int or not 1<=s['maxima']['passes']<=8:
                raise RuntimeError('saved directtarget exceeds frozen8cycle bound')
            counts=s.get('checks')
            required_counts={'nonzero_fail','missing_checkpoints','invalid_passes','primitive_below_floor','nonfinite',
                             'contraction','resolved_sign_reversal','seq_formula','round_formula',
                             'exceeded_physical','exceeded_normalized','exceeded_norm','exceeded_raw_excess','exceeded_background_excess'}
            if (not isinstance(counts,dict) or set(counts)!=required_counts or
                any(type(v) is not int or v<0 for v in counts.values())):
                raise RuntimeError('saved target reducer check counts missing/malformed')
            if any(counts.values()) or s['producer']['failed_trajectories']:
                raise RuntimeError('saved target check failure')
            if s['driver']['binary_sha256']!=binary_hash:raise RuntimeError('saved node uses different binary')
            if {b['sha256'] for b in s['source_bindings']}!=qualified_hashes:
                raise RuntimeError('saved target implementation sources differ')
            partials=s.get('partials')
            if (not isinstance(partials,list) or len(partials)!=9 or
                any(not isinstance(x,dict) for x in partials) or
                {(x['Nt'],x['Nw']) for x in partials}!={(nt,nw) for nt in (32,64,128) for nw in (128,256,512)}):
                raise RuntimeError('saved target nine quadratures incomplete or duplicated')
            meta_path=path.with_name(path.stem+'.execution.json')
            meta_rel=str(meta_path.relative_to(self.work))
            if not meta_path.is_file() or self.inputs.get(meta_rel)!=sha(meta_path):
                raise RuntimeError('missing/changed input-bound node CPU companion')
            meta=read(meta_path)
            if meta['identity']!=key or meta['summary_sha256']!=digest or meta['binary_sha256']!=binary_hash:
                raise RuntimeError('saved node CPU companion identity mismatch')
            for field in ('complete_local_cpu_seconds','complete_wall_seconds','component_self_plus_awaited_children_cpu_seconds','outer_controller_cpu_seconds_upper'):
                if field not in meta or type(meta[field]) not in (int,float) or not math.isfinite(meta[field]) or meta[field]<0:
                    raise RuntimeError('saved node complete CPU/wall companion field missing/malformed')
            if meta['complete_local_cpu_seconds']!=meta['component_self_plus_awaited_children_cpu_seconds']+meta['outer_controller_cpu_seconds_upper']:
                raise RuntimeError('saved directnode completeCPU components disagree')
            s['execution_measurements']=meta
            found[key]=s
        return found

    def forecast(self,records,done,rates):
        l=self.current_ledger();a=self.allowances()
        pending=[r for r in records if r['id'] not in done]
        main=sum(r['Nx']/256 for r in pending if r['family']=='main')
        diag=sum(r['Nx']/256 for r in pending if r['family']=='diagnostic')
        nm=sum(r['family']=='main' for r in pending);nd=len(pending)-nm
        C=float(l['charged_cpu_seconds'])+(cpu()-self.cpu_start)+outer_cpu_delta(outer_cpu_sample(self.ctx,self.active),self.outer_allocation_baseline)
        cp=C+1.5*main*rates['cpu']+max(308.3156323,1.5*(nm*a['per_node_host_cpu_seconds']+a['remaining_main_cpu']))
        cp+=max(1800.,1.5*(diag*rates['cpu']+nd*a['per_node_host_cpu_seconds']+a['remaining_diagnostic_cpu']))
        cp+=600+60+1.5*a['remaining_other_cpu']
        envelope=l['prospective_device_wall'];active_wall=time.monotonic()-float(self.active['allocation_start_monotonic'])
        wp=envelope['charged_science_wall_seconds']+active_wall+1.5*(main+diag)*rates['wall']
        wp+=1.5*(a['remaining_wall_seconds']+a['checkpoint_wall_seconds'])
        cert=read(self.asset('controller_certificate'))
        generated=float(cert['existing_generated_bytes'])+sum(x.stat().st_size for x in self.out.rglob('*') if x.is_file())
        # Count previous actual files and future outputs; historical baseline must be current bound, never reset.
        generated+=7864320+896*65536+float(cert['remaining_generated_bytes_allowance'])
        # All896node64KiB envelopes include prior/generated future summaries and
        # CPUcompanions. Current outputs may be conservatively counted again.
        if cp>14400 or wp>43200 or wp>57600 or generated>500000000:
            raise RuntimeError('complete remaining CPU/GPUupper/wall/artifact forecast fails')
        return dict(cpu=cp,whole_science_wall_gpu_upper=wp,artifact_bytes=generated,
                    remaining_main_equivalents=main,remaining_diagnostic_equivalents=diag,
                    remaining_nodes=len(pending),rates=rates.copy())

    def saved_controls(self):
        # All systematic inputs must come from bound SAVED sources, never a fresh primary scan.
        endpoint=read(self.asset('endpoint_certificate'))
        if not endpoint.get('pass') or endpoint['fine']['x']!=self.protocol['frozen_state']['x']:
            raise RuntimeError('saved source endpoint certificate mismatch/failure')
        if endpoint['fine']['amplitude_curvature']<=0 or endpoint['fine']['relaxed_quartic']<=0:
            raise RuntimeError('source endpoint amplitude/quartic control failed')
        wf=read(self.asset('frozen_witness'))
        for key in ('x','K','v','cP','EP','Eep'):
            if wf[key]!=self.protocol['frozen_state'][key]:raise RuntimeError('fixed witness binding differs')
        controls=read(self.asset('source_control_certificate'));Ectrl=finite(controls['Ectrl'])
        if not controls.get('pass_',controls.get('pass',False)):
            raise RuntimeError('bound source-control certificate did not pass')
        if Ectrl>float(self.protocol['frozen_state']['EP'])+float(self.protocol['frozen_state']['Eep']):
            raise RuntimeError('source control error exceeds inherited primary guard EP+Eep')
        primary=read(self.asset('primary_evaluations'));frozen=self.protocol['frozen_state']
        rows=[r for r in primary['records'] if r['x']==frozen['x'] and r['K']==frozen['K'] and '256,512' in r['table']]
        if not rows:raise RuntimeError('saved selected primary tail row missing')
        values=[]
        import numpy as np
        for row in rows:
            z=row['table']['256,512'];tail=np.array(z['tail5']['real'])+1j*np.array(z['tail5']['imag'])
            # Conservative existing operator-norm tail, plus saved primitive matrix error.
            values.append(2*float(np.linalg.norm(tail,2))+2*float(z['b_matrix']))
        shared_tail=max(values)
        if 'shared_tail_certificate' in self.assets:
            tc=read(self.asset('shared_tail_certificate'))
            if tc.get('fixed_K')!=frozen['K'] or tc.get('fixed_v')!=frozen['v'] or not tc.get('source_checks',{}).get('pass'):
                raise RuntimeError('saved shared-tail source/witness certificate mismatch')
            bound=tc.get('shared_tail_error',tc.get('error',tc.get('value')))
            if bound is None:raise RuntimeError('explicit saved shared-tail value missing')
            shared_tail=finite(bound)
            if shared_tail<max(values):
                raise RuntimeError('saved shared-tail certificate weaker than source-bound conservative norm')
        return shared_tail,Ectrl

    def target(self):
        frozen=self.protocol['frozen_state']
        self.guard();release=self.check_release();qual_call=int(release['qualification_call'])
        s=cached_state(self.root/'state.json');qc=s['calls'][qual_call-1]
        if qc['stage']!='science' or qc['version']!=self.ctx['version'] or qc['status']!='completed':
            raise RuntimeError('same-version successful actual qualification call required')
        qpath=self.prior_file(qual_call,release['qualification_result_sha256']);qualified_record=read(qpath)
        qr=qualified_record.get('report',qualified_record)
        if qr.get('candidate_sha256')!=sha(self.candidate_path):raise RuntimeError('qualification report candidate differs')
        if qr.get('qualification_pass') is not True or qr.get('target_nodes',0)!=0:
            raise RuntimeError('qualification was not a target-free pass')
        self.prior_file(qual_call,release['qualification_manifest_sha256'])
        source=self.prior_file(qual_call,release['qualified_source_sha256'])
        if sha(source)!=self.assets['direct_cuda_source']['sha256'] or source.read_bytes()!=self.asset('direct_cuda_source').read_bytes():
            raise RuntimeError('qualified directsource differs from exact admitted source')
        binary=self.prior_file(qual_call,release['qualified_binary_sha256'])
        rates=dict(cpu=finite(qr['measurements']['stress_complete_local_cpu_seconds']),
                   wall=finite(qr['measurements']['stress_wall_seconds']),
                   kernel=finite(qr['measurements']['kernel_event_seconds']))
        if min(rates.values())<=0:raise RuntimeError('qualified rate missing')
        # Mandatory source/endpoint/frozen-witness controls before ANY target node, not only final assembly.
        shared_tail,Ectrl=self.saved_controls()
        driver=load_module(self.asset('driver'),'target_driver');reducer=load_module(self.asset('reducer'),'target_reducer')
        reduction=driver.original_reduction(self.asset('sums_library'));records=driver.schedule()
        if len(records)!=896:raise RuntimeError('frozen full schedule count changed')
        expected={r['id']:r for r in records}
        done=self.validated_prior_nodes(expected,sha(source),sha(binary))
        # Never decrease qualified or observed worst rates for later faster nodes.
        for r in records:
            if r['id'] in done:
                v=done[r['id']];scale=r['Nx']/256
                rates['wall']=max(rates['wall'],finite(v['execution_measurements']['complete_wall_seconds'])/scale)
                rates['kernel']=max(rates['kernel'],v['producer']['kernel_seconds']/scale)
                measured=v.get('execution_measurements',{}).get('complete_local_cpu_seconds')
                if measured is None:raise RuntimeError('saved node complete CPU rate missing')
                rates['cpu']=max(rates['cpu'],float(measured)/scale)
        forecast=self.forecast(records,done,rates)
        source_bindings=[dict(path=str(x),sha256=sha(x)) for x in
                         (source,binary,self.asset('driver'),self.asset('reducer'),self.asset('sums_library'))]
        # Group immutable predeclared schedule into complete signed-node pairs.
        pairs=[records[i:i+2] for i in range(0,len(records),2)];used=0.;stopped=False
        folder=self.out/'nodes';folder.mkdir(exist_ok=True)
        for pair in pairs:
            if len(pair)!=2 or pair[0]['signed_node_id']!=pair[1]['signed_node_id']:
                raise RuntimeError('invalid frozen signed-node pair')
            flags=[r['id'] in done for r in pair]
            if any(flags) and not all(flags):raise RuntimeError('partial prior signed-node pair requires independent disposition, no automatic replay')
            if all(flags):continue
            eq=2*pair[0]['Nx']/256
            if used+eq>54:stopped=True;break
            estimated_artifact_rows=sum(x.is_file() for x in self.out.rglob('*'))+8
            if estimated_artifact_rows*220>60000:stopped=True;break
            self.guard();forecast=self.forecast(records,done,rates)
            if shutil.disk_usage(self.work).free<256*1024*1024:
                raise RuntimeError('disk free below256MiB; exact checkpoint retained, no archive/input deletion')
            needed=1.5*rates['wall']*eq+5
            if time.monotonic()+needed>self.deadline-90:stopped=True;break
            pair_new=[]
            for r in pair:
                self.guard();before=cpu();outer_before=outer_cpu_sample(self.ctx,self.active);then=time.monotonic()
                summary=driver.run_node(binary,r,self.protocol['frozen_state'],folder,reducer,reduction,
                    controller_permit=self.guard,timeout=min(3600.,self.deadline-time.monotonic()-45),
                    source_bindings=source_bindings,finalization_reserve=45)
                component_cpu=cpu()-before;outer_part=outer_cpu_delta(outer_cpu_sample(self.ctx,self.active),outer_before);dt=component_cpu+outer_part;wall=time.monotonic()-then
                if summary['driver']['parent_cpu_seconds']>self.allowances()['per_node_host_cpu_seconds']:
                    raise RuntimeError('observed target host/reducer CPU exceeds reviewed1.2s allowance')
                # Append metering in a separate bounded immutable companion; do
                # not edit the driver summary already committed with its hash.
                meta=dict(identity=r['id'],summary_sha256=sha(folder/(r['id']+'.json')),
                          complete_local_cpu_seconds=dt,component_self_plus_awaited_children_cpu_seconds=component_cpu,outer_controller_cpu_seconds_upper=outer_part,complete_wall_seconds=wall,
                          source_sha256=sha(source),binary_sha256=sha(binary))
                new_json(folder/(r['id']+'.execution.json'),meta)
                if (folder/(r['id']+'.execution.json')).stat().st_size+(folder/(r['id']+'.json')).stat().st_size>65536:
                    raise RuntimeError('combined node summary/CPUcompanion exceeds64KiB')
                summary['execution_measurements']=meta;pair_new.append((r,summary))
                scale=r['Nx']/256;rates['cpu']=max(rates['cpu'],dt/scale)
                rates['wall']=max(rates['wall'],wall/scale)
                rates['kernel']=max(rates['kernel'],summary['producer']['kernel_seconds']/scale)
                used+=scale
            for r,summary in pair_new:done[r['id']]=summary;self.new_summaries.append(summary)
            self.guard();forecast=self.forecast(records,done,rates)
        if len(done)==896:
            tables=driver.assemble_node_tables([done[r['id']] for r in records],frozen,reduction,driver.analytic_factory(frozen,reduction))
            certificate=driver.independent_error_assembly(tables,frozen,shared_tail,Ectrl)
            new_json(self.out/'INDEPENDENT_CERTIFICATE.json',certificate)
            new_json(self.out/'ASSEMBLED_TABLES.json',[dict(key=list(k),**v) for k,v in tables.items()])
            accepted=bool(certificate['numerically_accepted']);self.terminal=True
            return dict(mode='target',outcome='supported_empirical_fixed_witness' if accepted else 'inconclusive',
                        completed_target_nodes=896,new_target_nodes=len(self.new_summaries),
                        forecast=forecast,independent_certificate=certificate,
                        next_step='Actual independent result/parent coverage review; empirical envelope only, no automatic scientific acceptance or parent closure.')
        return dict(mode='target',outcome='inconclusive_chunk',completed_target_nodes=len(done),
                    new_target_nodes=len(self.new_summaries),forecast=forecast,
                    next_step='Independent result review of this bounded chunk, then root continuation only for missing complete signed-node pairs; no replay or rate reduction.')

    def finish(self,report,error=None):
        if error:report=dict(mode=self.args.mode,outcome='inconclusive',blockers=[type(error).__name__+': '+str(error)],
                             new_target_nodes=len(list((self.out/'nodes').glob('*.json')))//2 if (self.out/'nodes').exists() else 0,
                             next_step='Independent result review/controller disposition; partial signed-node failure cannot automatically resume or replay.')
        report.update(disk_free_bytes_observed=shutil.disk_usage(self.work).free,
                      disk_scope='Per-pair256MiB guard supplements root pre-claim immutable-copy projection; no completeness guarantee and no deletion authorized.',
                      candidate_sha256=sha(self.candidate_path),execution_context_sha256=sha(self.work/'execution_context.json'),
                      complete_adapter_local_cpu_seconds=cpu()-self.cpu_start,adapter_wall_seconds=time.monotonic()-self.start,
                      local_self_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      awaited_child_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                      measurement_limitations='Adapter SELF+awaited CHILD CPU excludes outer controller/guard/input/intake; root adds actual outer controller SELF+CHILD. Separate unsynchronized RSS maxima; root sampled process-group and whole-sciencewall GPUupper accounting authoritative. GPUevent rate is not complete GPUcontext memory or consumed CPU. Exit/receipt costs remain charged by root.',
                      no_case_or_ledger_mutations=True)
        new_json(self.out/'RESULT.json',report)
        artifacts=[str(x.relative_to(self.work)) for x in sorted(self.out.rglob('*')) if x.is_file()]
        new_json(self.out/'MANIFEST.json',dict(artifacts=[dict(path=x,sha256=sha(self.work/x),bytes=(self.work/x).stat().st_size) for x in artifacts]))
        artifacts.append('science/MANIFEST.json')
        result=dict(binding_sha256=self.ctx['binding_sha256'],decision='result',summary='Deterministic '+self.args.mode+': '+report['outcome']+'. No automatic target release or parent closure.',
                    findings=['Actual hash-bound external deterministic '+self.args.mode+' computation; see saved science/RESULT.json.',
                              'All historical failures/timeouts/charges preserved; no inspected primary data relabeled fresh.',
                              'Outcome '+report['outcome']+' requires independent result review.'],
                    blockers=report.get('blockers',[]) or report.get('report',{}).get('blockers',[]),
                    local_exhaustion=[],next_step=report['next_step'],artifacts=artifacts,
                    progress_evidence=[dict(change='Saved actual '+self.args.mode+' evidence and bounded accounting.',reference='science/RESULT.json')],
                    candidate_json=None,novelty_review_json=None,external_dependency=None,question_assessment_json=None)
        Path(self.args.result).write_text(json.dumps(result,indent=2)+'\n')
        return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',required=True)
    p.add_argument('--case',required=True);p.add_argument('--mode',choices=('qualification','target'),required=True)
    p.add_argument('--result',required=True);args=p.parse_args()
    adapter=None
    try:
        adapter=Adapter(args)
        report=adapter.qualification() if args.mode=='qualification' else adapter.target()
        adapter.finish(report);return 0
    except BaseException as error:
        if adapter is not None:
            adapter.finish(None,error);return 0 # genuine valid inconclusive result, not technical exit
        # Gate/input preparation failure before any science remains technical.
        print(type(error).__name__+': '+str(error),file=sys.stderr);return 2


if __name__=='__main__':raise SystemExit(main())
