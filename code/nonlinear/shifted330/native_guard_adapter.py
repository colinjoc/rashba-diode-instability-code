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
        measured_cpu = cpu()-self.cpu_start+outer_cpu_delta(outer_cpu_sample(self.ctx,self.active),self.outer_allocation_baseline)
        cap = l['active_allocation']['per_invocation_cpu_cap_seconds']
        policy = load_module(self.asset('cpu_policy'), 'native_exact_cpu_policy')
        effective = policy.runtime_limit(self.root,self.candidate,read(self.root/'orchestrator/SCIENCE_RELEASE.json'),l)
        tasks = policy.finite_task_limits(self.protocol)
        if cap is None:
            if effective is not None or tasks['invocation_cpu_seconds'] is not None or tasks['qualification_cpu_seconds'] is not None or not self.protocol.get('cpu_cap_removal_authorization'):
                raise RuntimeError('Null invocation CPU requires exact user/admission/release policy')
        elif policy.cpu_exceeded(measured_cpu,cap):
            raise RuntimeError('local CPU allowance exhausted')
        for name in self.assets:self.asset(name)
        if self.args.mode=='target':self.check_release()

    def check_release(self):
        r=read(self.root/'orchestrator/TARGET_EXECUTION_RELEASE.json')
        if r.get('scope')!='prospective_shifted_angular_fixed_witness_confirmation' or r.get('version')!=self.ctx['version'] or r.get('call_id')!=self.ctx['call_id'] or r.get('candidate_sha256')!=sha(self.candidate_path):
            raise RuntimeError('exact target release missing/changed')
        if r.get('protocol_sha256')!=sha(self.asset('qualification_protocol')) or r.get('worker_sha256')!=sha(self.asset('native_worker')):
            raise RuntimeError('exact target protocol/worker release changed')
        cc=read(self.work/'confirmation_context.json')
        confirmation=self.protocol['confirmation']
        identities=confirmation.get('identities', [])
        actual=cc.get('this_call')
        if (confirmation.get('mode')!='heldout' or len(identities)!=1 or
                r.get('confirmation_identity')!=identities[0] or
                not isinstance(actual,dict) or actual.get('mode')!='heldout' or
                actual.get('reserved')!=identities[0]):
            raise RuntimeError('exact heldout target confirmation identity/release/reservation required')
        return r








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


