c = get_config()
c.Device.zmq_port=12345
c.Modeler.use_custom_metric = True
c.Modeler.custom_metric_path = 'metric.pickl'
c.Modeler.clusterer = 'hybrid'
c.OpenMMSimulator.report_interval = 10
c.OpenMMSimulator.number_of_steps = 2400000
c.OpenMMSimulator.minimize = False
c.AdaptiveServer.topology_pdb = '../../../single.pdb'
c.BaseSampler.seed_structures = 'seed_structures.h5'
c.CountsSampler.beta = 1
c.Modeler.lag_time = 10
