from pathlib import Path
import numpy as np
import pytest
from pyruvseq import ruvr, ruvs, make_groups

DATA=Path(__file__).parents[1]/'data'
@pytest.mark.parametrize('case',['r','r_log','s','s_log'])
def test_r_reference(case):
    x=np.loadtxt(DATA/'input.csv',delimiter=',')
    is_log=case.endswith('log')
    x=np.log1p(x) if is_log else x
    options=dict(round_counts=False,is_log=is_log)
    if case.startswith('r'):
        out=ruvr(x,range(12),2,np.loadtxt(DATA/'residuals.csv',delimiter=','),center=not is_log,**options)
    else:
        groups=make_groups(np.repeat(list('abcd'),4))
        np.testing.assert_array_equal(groups+1,np.loadtxt(DATA/'groups_r.csv',delimiter=','))
        out=ruvs(x,range(12),2,groups,**options)
    w=out['W']; projection=w@np.linalg.solve(w.T@w,w.T)
    np.testing.assert_allclose(projection,np.loadtxt(DATA/f'{case}_projection.csv',delimiter=','),atol=1e-8,rtol=1e-8)
    np.testing.assert_allclose(out['normalized_counts'],np.loadtxt(DATA/f'{case}_corrected.csv',delimiter=','),atol=1e-8,rtol=1e-8)

def test_invalid_replicates_and_residuals():
    x=np.loadtxt(DATA/'input.csv',delimiter=',')
    with pytest.raises(ValueError,match='same shape'): ruvr(x,[0,1],1,np.ones((2,2)))
    with pytest.raises(ValueError,match='two samples'): ruvs(x,[0,1],1,np.arange(16)[:,None])
    with pytest.raises(ValueError,match='repeat'): ruvs(x,[0,1],1,[[0,0]])
