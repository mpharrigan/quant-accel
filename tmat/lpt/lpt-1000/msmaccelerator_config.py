c = get_config()
c.Modeler.clusterer = 'kcenters'
c.Modeler.rmsd_atom_indices = '../../AllAtoms.dat'
c.AdaptiveServer.md_engine = 'TMat'
c.AdaptiveServer.topology_pdb = '../../ntl9_heavy.pdb'
c.AdaptiveServer.gens_fn = '../../Gens.h5'
c.BaseSampler.seed_structures = '../../Gens.h5'

# TMat stuff
c.TMatSimulator.tmat_fn = '../../tProb.mtx'
c.TMatSimulator.gens_fn = '../../Gens.h5'
c.TMatSimulator.report_interval = 1

# Run dependent stuff
c.TMatSimulator.number_of_steps = 1000
c.Modeler.lag_time = 50
c.Modeler.stride = 1
c.BaseServer.zmq_port = 12342
c.Device.zmq_port = 12342
