c = get_config()
c.Device.zmq_port = 12349
c.Modeler.clusterer = 'kcenters'
c.Modeler.lag_time = 50
c.Modeler.rmsd_atom_indices = '../../AllAtoms.dat'
c.BaseServer.zmq_port = 12349
c.AdaptiveServer.md_engine = 'TMat'
c.AdaptiveServer.topology_pdb = '../../ntl9_heavy.pdb'
c.AdaptiveServer.gens_fn = '../../Gens.h5'
c.BaseSampler.seed_structures = '../../Gens.h5'

# TMat stuff
c.TMatSimulator.tmat_fn = '../../tProb.mtx'
c.TMatSimulator.gens_fn = '../../Gens.h5'
c.TMatSimulator.report_interval = 1

# Run dependent stuff
c.TMatSimulator.number_of_steps = 100
