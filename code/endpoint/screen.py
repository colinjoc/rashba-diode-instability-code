import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,resource
from pathlib import Path
import numpy as np
from primary import Primary,certificate,decode,encode,scalar_table
out=Path('science');ep=json.loads((out/'endpoint_certificate.json').read_text());br=json.loads((out/'branch_topology.json').read_text())
assert ep['pass'] and br['pass']
x=np.array(ep['fine']['x']);T,D,q,h=x
p=Primary();rows=[];eligible=[]
try:
    for lam in [-1,1]:
        k0=2*lam*h-q
        for angle in [0,np.pi/4,np.pi/2,3*np.pi/4]:
            for scale in [.5,.75,1,1.25,1.5]:
                K=k0*scale*np.array([np.cos(angle),np.sin(angle)])
                if not .25<=np.linalg.norm(K)<=20:
                    rows.append({'K':K.tolist(),'excluded':'outside frozen prepared-code range'});continue
                # Three prescribed screenings; the fine evaluation retains its nested integrands.
                a=p.evaluate(x,K,64,128,nested=False)
                b=p.evaluate(x,K,128,256,nested=False)
                c=p.evaluate(x,K)
                cert=certificate(c)
                m=decode(c['table']['256,512']['matrix']);vals,vecs=np.linalg.eigh(m);v=vecs[:,0]
                center=float(vals[0]);err=cert['EP']
                row={'lambda':lam,'angle':angle,'scale':scale,'K':K.tolist(),'center':center,'eigenvalues':vals.tolist(),
                     'v':encode(v),'EP':err,'certificate':cert,'direction_curvatures':scalar_table(c,v),
                     'nested_direct_discrepancy':max(float(np.max(abs(decode(a['table']['64,128']['matrix'])-decode(c['table']['64,128']['matrix'])))),
                         float(np.max(abs(decode(b['table']['128,256']['matrix'])-decode(c['table']['128,256']['matrix'])))))}
                row['eligible_primary']=bool(cert['pass'] and center+err<0 and err<=.25*abs(center))
                rows.append(row)
                if row['eligible_primary']:eligible.append(row)
                print('MODE',len(rows),K.tolist(),center,err,cert['pass'],row['eligible_primary'],flush=True)
    # Keep selection prospective to nonlinear data. Tie handling uses overlapping primary error envelopes.
    chosen=None
    if eligible:
        best=min(eligible,key=lambda r:r['center'])
        ties=[r for r in eligible if abs(r['center']-best['center'])<=r['EP']+best['EP']]
        chosen=min(ties,key=lambda r:tuple(r['K']))
    result={'status':'primary exploratory screening only; no nonlinear confirmation','rows':rows,'selected_provisional':chosen,
            'matrix_evaluations':p.count,'primary_cpu_seconds':time.process_time()-p.start,
            'local_process_maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'All 40 frozen candidate vectors; no local search expansion or K-to-zero witness.'}
    (out/'primary_screen.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
finally:
    (out/'primary_evaluations.json').write_text(json.dumps({'records':p.records,'count':p.count,'cpu_seconds':time.process_time()-p.start},indent=2,allow_nan=False)+'\n')
