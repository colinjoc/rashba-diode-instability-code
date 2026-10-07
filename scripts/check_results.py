"""Audit saved scientific margins; never generate fresh solver measurements."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def check():
    data = json.loads((ROOT / 'scientific_data.json').read_text())
    linear = json.loads((ROOT / 'data/evidence/LINEAR271.json').read_text())
    nonlinear = json.loads((ROOT / 'data/evidence/qualification330.json').read_text())
    f = data['perturbation']
    cp, ep, eep = f['cP'], f['EP'], f['Eep']
    ci, ei = linear['independent_coefficient'], linear['complete_physical_EI']
    checks = {
        'saved_linear_pass': linear['qualification_pass'] and not linear['failed_checks'],
        'saved_nonlinear_pass': nonlinear['whole_qualification_pass'] and not nonlinear['failed_checks'],
        'primary_negative': cp + ep + eep < 0,
        'independent_negative': ci + ei + eep < 0,
        'primary_quarter_error': ep + eep <= .25 * abs(cp),
        'independent_quarter_error': ei + eep <= .25 * abs(ci),
        'primary_independent_agreement': abs(cp - ci) <= ep + ei,
        'control_adequacy': linear['Ectrl'] <= min(ep + eep, ei + eep),
        'published_linear_matches_saved': data['linear_response']['independent_coefficient'] == ci and data['linear_response']['complete_physical_EI'] == ei,
        'published_nonlinear_matches_saved': data['shifted_angular']['complete_results'] == nonlinear['complete_results'] and data['shifted_angular']['rows'] == nonlinear['rows'],
        'complete_schedule': len(nonlinear['rows']) == 93 and len({tuple(r['configuration']) for r in nonlinear['rows']}) == 31,
        'angular_origins_agree': all(r['pass_'] and r['difference'] <= r['allowed'] for r in nonlinear['origin_comparisons']),
    }
    rows = []
    for index, row in enumerate(nonlinear['complete_results']):
        upper = row['coefficient'] + row['EI'] + eep
        checks[f'finite_amplitude_{index}_negative'] = upper < 0
        rows.append(dict(epsilon=row['epsilon'], coefficient=row['coefficient'], numerical_allowance=row['EI'], endpoint_allowance=eep, upper_margin=upper))
    report = dict(scope='saved-evidence audit; no new confirmation', primary_coefficient=cp,
                  independent_coefficient=ci, independent_numerical_allowance=ei,
                  endpoint_allowance=eep, nonlinear=rows, checks=checks)
    print(json.dumps(report, indent=2))
    if not all(checks.values()):
        raise RuntimeError('Saved evidence audit failed: ' + ', '.join(k for k, v in checks.items() if not v))
    return report

if __name__ == '__main__':
    check()
