import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,itertools,hashlib,resource
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from uniform import Uniform,eseq
from primary import Primary,certificate,decode,encode,norm
out=Path('science')
ep=json.loads((out/'endpoint_certificate.json').read_text());sc=json.loads((out/'primary_screen.json').read_text())
old=json.loads((out/'primary_evaluations.json').read_text());p=Primary();p.count=old['count'];p.records=old['records'];p.start=time.process_time()-old['cpu_seconds']
ud=json.loads((out/'uniform_evaluations.json').read_text());u=Uniform(cpu_cap=500);u.count=ud['count'];u.log=ud['evaluations'];u.start=time.process_time()-ud['cpu_seconds']
x=np.array(ep['fine']['x']);rho=np.array(ep['rho']);sel=sc['selected_provisional'];assert sel
K=np.array(sel['K']);v=decode(sel['v']);T,D,q,h=x
details={}
try:
    anchor=p.evaluate(x,K,cost=9);cp=certificate(anchor)
    center=float(np.vdot(v,decode(anchor['table']['256,512']['matrix'])@v).real)
    def scalar(xx):
        r=p.evaluate(xx,K);c=certificate(r)
        return {'x':list(xx),'c':float(np.vdot(v,decode(r['table']['256,512']['matrix'])@v).real),'EP':c['EP'],'pass':c['pass']}
    axes=[];s=[];sb=[];g=[];gb=[]
    for i in range(4):
        probes=[scalar(x+a*rho[i]*np.eye(4)[i]) for a in [-1,-.5,.5,1]]
        vals=np.array([a['c'] for a in probes]);errs=np.array([a['EP'] for a in probes])
        slopes=(vals-center)/(np.array([-1,-.5,.5,1])*rho[i]);s.append(max(abs(slopes)))
        bounds=(errs+cp['EP'])/(np.array([1,.5,.5,1])*rho[i]);sb.append(max(bounds))
        full=(vals[3]-vals[0])/(2*rho[i]);half=(vals[2]-vals[1])/rho[i]
        bf=(errs[3]+errs[0])/(2*rho[i]);bh=(errs[2]+errs[1])/rho[i]
        ok=abs(full-half)<=.25*max(abs(full),abs(half))+bf+bh
        g.append(half);gb.append(bh)
        axes.append({'i':i,'probes':probes,'one_sided_slopes':slopes.tolist(),'slope_bounds':bounds.tolist(),
                     'full_centered':full,'half_centered':half,'smoothness_pass':bool(ok)})
    corners=[];remainders=[];rembounds=[]
    for signs in itertools.product([-1,1],repeat=4):
        dx=rho*np.array(signs);r=scalar(x+dx)
        rem=r['c']-center-np.dot(g,dx);bound=r['EP']+cp['EP']+np.dot(gb,abs(dx))
        r.update(signs=signs,linear_remainder=rem,remainder_bound=bound);corners.append(r)
        remainders.append(abs(rem));rembounds.append(bound)
    enom=2*np.dot(rho,s)+2*max(remainders)
    eprop=2*np.dot(rho,sb)+2*max(rembounds)
    eep=enom+eprop
    epok=all(a['smoothness_pass'] for a in axes) and all(a['pass'] for a in corners) and all(b['pass'] for a in axes for b in a['probes'])
    details['endpoint_propagation']={'pass':bool(epok),'axes':axes,'corners':corners,'Eep_nominal':enom,'Eep_error_propagation':eprop,'Eep':eep}
    # Source Ward and zero momentum amplitude control (raw matrix).
    zero=p.evaluate(x,[0,0]);m0=decode(zero['table']['256,512']['matrix']);e0=certificate(zero)['EP']
    uniform=u.evaluate(x,256,512);fd=uniform['F_derivatives']
    ward_expected=2*uniform['g'][0]/D
    warderr=abs(m0[1,1].real-ward_expected)
    ampdefect=abs(m0[0,0].real-uniform['amplitude_curvature'])
    wardok=warderr<=1e-8 and ampdefect<=1e-8 and m0[0,0].real-e0>0
    # Small-K identities with prescribed Richardson rules and propagated errors.
    ks=[.01,.005,.0025];cross=[];phase=[];bc=[];bp=[];small=[]
    for k in ks:
        rr=p.evaluate(x,[k,0]);ct=certificate(rr);m=decode(rr['table']['256,512']['matrix']);e=ct['EP']
        c=m[0,1].imag/k;schur=m[1,1].real-abs(m[0,1])**2/m[0,0].real
        ph=D*D*schur/(k*k)
        bphase=D*D/(k*k)*(e+2*abs(m[0,1])*e/(m[0,0].real-e)+abs(m[0,1])**2*e/((m[0,0].real-e)*m[0,0].real)+e*e/(m[0,0].real-e))
        cross.append(c);bc.append(e/k);phase.append(ph);bp.append(bphase)
        small.append({'K':k,'mixed':c,'phase_Schur':ph,'mixed_error':e/k,'phase_error':bphase,'kernel_pass':ct['pass']})
    def rich(y,b,expected):
        err,ok=eseq(y,b);r01=(4*y[1]-y[0])/3;r12=(4*y[2]-y[1])/3
        br01=(4*b[1]+b[0])/3;br12=(4*b[2]+b[1])/3
        extra=4*(abs(r12-r01)+br01+br12)+br12
        return {'pass':bool(ok and abs(r12-expected)<=1e-8+extra),'raw_sequence_pass':bool(ok),
                'expected':expected,'R01':r01,'R12':r12,'bound':extra,'difference':abs(r12-expected),'Eseq':float(err)}
    mix=rich(cross,bc,fd['1,1']/D);schur=rich(phase,bp,uniform['g'][2]/2)
    # Zero-field BCS recovery and equality to prepared source functional.
    def gap(d):return u.evaluate([.5,d,0,0],128,256)['g'][0]
    db=brentq(gap,1,2,xtol=1e-13);ub=u.evaluate([.5,db,0,0],256,512)
    bcsmat=p.evaluate([.5,db,0,0],[0,0]);bcsm=decode(bcsmat['table']['256,512']['matrix'])
    bcs={'T':.5,'Delta':db,'gap_residual':ub['g'][0],'current':ub['g'][1],'phase_kernel':bcsm[1,1].real,
         'amplitude_kernel':bcsm[0,0].real,'pass':bool(abs(ub['g'][0])<1e-8 and abs(ub['g'][1])<1e-8 and abs(bcsm[1,1])<1e-8 and bcsm[0,0].real>0)}
    controls={'ward_error':warderr,'ward_expected':ward_expected,'K0_phase':m0[1,1].real,
         'amplitude_difference':ampdefect,'amplitude_kernel':m0[0,0].real,'amplitude_error':e0,
         'ward_and_amplitude_pass':bool(wardok),'smallK':small,'mixed_identity':mix,'schur_identity':schur,'zero_field_BCS':bcs,
         'raw_hermiticity_max':max(t['hermiticity_defect'] for r in p.records for t in r['table'].values())}
    controls['pass']=bool(wardok and mix['pass'] and schur['pass'] and bcs['pass'] and controls['raw_hermiticity_max']<=1e-12)
    valid=bool(epok and cp['pass'] and controls['pass'] and center+cp['EP']+eep<0 and cp['EP']+eep<=.25*abs(center))
    details.update(status='primary selected-witness checks',primary_pass=valid,controls=controls,
         x=x.tolist(),K=K.tolist(),v=encode(v),cP=center,EP=cp['EP'],Eep=eep,primary_upper=center+cp['EP']+eep,
         relative_error=(cp['EP']+eep)/abs(center),matrix_evaluations=p.count,
         uniform_evaluations=u.count,primary_cpu_seconds=time.process_time()-p.start,uniform_cpu_seconds=time.process_time()-u.start)
    (out/'selected_checks.json').write_text(json.dumps(details,indent=2,allow_nan=False)+'\n')
    if valid:
        frozen={k:details[k] for k in ['x','K','v','cP','EP','Eep','primary_upper']}
        frozen.update(test_id='rashba-source-endpoint-fixed-finiteK-paired-curvature-v2',
          candidate_sha256='12621a3ea333099972beae5b6b2fa370c0b842edbba1f85d46887826cae17d25',
          status='Fixed before any nonlinear perturbation data; no replacement allowed',
          epsilon=[.01,.005,.0025],created_epoch=time.time(),created_monotonic=time.monotonic())
        (out/'frozen_witness.json').write_text(json.dumps(frozen,indent=2,allow_nan=False)+'\n')
    print('SELECTED',json.dumps({k:v for k,v in details.items() if k not in ['controls','endpoint_propagation']}),flush=True)
    print('CONTROLS',json.dumps(controls),flush=True)
finally:
    (out/'primary_evaluations.json').write_text(json.dumps({'records':p.records,'count':p.count,'cpu_seconds':time.process_time()-p.start},indent=2,allow_nan=False)+'\n')
    u.save(out/'uniform_evaluations.json',{'phase':'selected-witness controls'})
