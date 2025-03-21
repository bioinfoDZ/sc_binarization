import muon as mu


from muon import atac as ac
import cupy as cp
import omicverse as ov
import scanpy as sc


rna = sc.read('sim_rna.h5ad')
atac =sc.read('sim_atac.h5ad')

test_mofa = ov.single.pyMOFA(omics=[rna, atac], omics_name=["RNA", "ATAC"])

test_mofa.mofa_preprocess()

test_mofa.mofa_run(gpu_mode=True, outfile='sim_mofa.hdf5')
