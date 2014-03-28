'''
Created on Mar 5, 2014

@author: harrigan

Generate system and integrator files for running openmm.
This can be run from the command line
'''

from simtk import openmm, unit
import toy_accel.mullerforce as mf
import argparse
import logging as log


def generate_openmm(sys_fn, int_fn):
    mass = 12.0 * unit.dalton
    temperature = 750 * unit.kelvin
    friction = 100 / unit.picosecond
    timestep = 10.0 * unit.femtosecond

    # Prepare the system
    system = openmm.System()
    mullerforce = mf.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    log.info("Writing System file: %s", sys_fn)
    with open(sys_fn, 'w') as f:
        f.write(openmm.XmlSerializer.serialize(system))

    # Prepare integrator
    log.info("Writing Integrator file: %s", int_fn)
    integrator = openmm.LangevinIntegrator(temperature, friction, timestep)
    with open(int_fn, 'w') as f:
        f.write(openmm.XmlSerializer.serialize(integrator))


def parse():
    parser = argparse.ArgumentParser(description='Perform accelerated sampling',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sys_fn', help='System file', default='system.xml')
    parser.add_argument(
        '--int_fn',
        help='Integrator file',
        default='integrator.xml')
    args = parser.parse_args()
    generate_openmm(args.sys_fn, args.int_fn)


if __name__ == "__main__":
    log.basicConfig(level=log.INFO)
    parse()
