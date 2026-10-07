"""Exact-import synthetic algebra/derivative fixtures only, never physical endpoint."""
import importlib.util,json,math,sys,time,resource
from pathlib import Path
B=Path(__file__).resolve().parent;sys.dont_write_bytecode=True
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
m=load(B/'integrated_source_budget.py','synthetic_integrated');lin=load(B/'linear_response.py','synthetic_linear')
def point(p,vals):
 return sum(v*math.prod(vals[j]**e[j] for j in range(4)) for e,v in p.items())
start=time.process_time();wall=time.monotonic();rows=[]
v=(complex(.3,.1),complex(-.2,.5));K=.4;W,Q,N,polys=m.formulas(K,v)
for vals in [(1.1,.8,.6,.2),(.3,1.4,-.7,.85),(3.,.4,.5,-.4)]:
 w,D,a,c=vals;z=complex(w,a*c);omega=(z*z+D*D)**.5;den=omega*(4*(z*z+D*D)+(K*c)**2)
 contracted=(point(N,vals)/den).real;_,H=lin.force_matrix(z,D,K*c);reference=lin.contract(H,v);err=abs(contracted-reference)
 derivatives=[]
 for index in (0,1,2):
  step=1e-6;lo=list(vals);hi=list(vals);lo[index]-=step;hi[index]+=step
  def S(x):
   z=complex(x[0],x[2]*x[3]);return (point(N,x)/((z*z+x[1]*x[1])**.5*(4*(z*z+x[1]*x[1])+(K*x[3])**2))).real
  fd=(S(hi)-S(lo))/(2*step)
  expected=contracted+w*fd if index==0 else fd
  box=tuple((x-1e-7,x+1e-7) for x in vals);enclosure,_,_=m.real_ratios(polys,box,K);interval=enclosure[index]
  derivatives.append({'variable':['joint_T','D','a'][index],'finite_difference':expected,'enclosure':list(interval),'pass':interval[0]-1e-8<=expected<=interval[1]+1e-8})
 rows.append({'synthetic_values':vals,'contract_discrepancy':err,'derivatives':derivatives,'pass':err<1e-12 and all(x['pass'] for x in derivatives)})
eb=load(B/'error_bounds.py','synthetic_tail')
synthetic_witness={'K':[K,0.],'v':{'real':[x.real for x in v],'imag':[x.imag for x in v]}}
synthetic_source=[.7,1.2,.3,.1];synthetic_radii=[1e-7]*4
tail=eb.corrected_tail(synthetic_source,8,synthetic_witness,synthetic_radii)
synthetic_full=m.continuous_source_budget(synthetic_source,synthetic_radii,8,16,synthetic_witness,tail,lambda:None)
report={'synthetic_complete_budget':synthetic_full,'scope':'Synthetic exact imported contracted algebra/real denominator/derivative production path; no endpoint bounds or physical responses','physical_execution':False,'endpoint_evaluations':0,'all_pass':all(x['pass'] for x in rows),'rows':rows,'author_fixture_self_cpu_seconds':time.process_time()-start,'author_fixture_wall_seconds':time.monotonic()-wall,'lifetime_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
print(json.dumps(report,sort_keys=True,allow_nan=False))
if not report['all_pass']:raise SystemExit(1)
