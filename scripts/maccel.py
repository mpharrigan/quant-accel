__author__ = 'harrigan'

import argparse
import logging
import os

import maccelerator as maccel


logging.basicConfig(level=logging.INFO)


def config_entry(args):
    """Entry point for argparse."""
    print('Writing some sample configuration files. Make sure you look at')
    print('them before blindly executing! They will need modification.')
    config(args.out_fn, args.n_copy, args.cluster, args.config, args.parallel)


def config(out_fn, n_copy, cluster, config_name, parallel):
    cluster_dict = dict(pbs=maccel.PBSCluster)
    config_dict = dict(alanine=maccel.AlanineConfiguration)
    grid_manager_name = 'MaccelGridShm'  # TODO: Make general

    clust = cluster_dict[cluster](n_copy=n_copy, parallel=parallel)
    config = config_dict[config_name]

    job_out = "{}.{}".format(out_fn, clust.job_script_ext)
    py_out = "{}.py".format(out_fn)

    with open(job_out, 'w') as job_f:
        job_f.write(clust.make_job_script(py_out))

    with open(py_out, 'w') as py_f:
        py_f.write(config.get_template(grid_manager_name))


def plot_entry(args):
    plot(args.run_fn, parallel=False)


def plot(run_fn, parallel):
    load_dir = os.path.dirname(run_fn)
    run = maccel.special_pickle_load(run_fn)
    pm = maccel.PlotMaker(run, load_dir=load_dir, parallel=parallel)
    pm.make_plots()


def parse():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sp = parser.add_subparsers(dest='command')
    sp.required = True

    plot_p = sp.add_parser('plot',
                           formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    plot_p.set_defaults(func=plot_entry)
    plot_p.add_argument('--run_fn', '-r',
                        help="""Path to run.pickl to make plots for.""",
                        default="run.pickl")



    # ----------------------------------------------------
    config_p = sp.add_parser('config',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    config_p.set_defaults(func=config_entry)
    config_p.add_argument('--cluster', '-c',
                          help="""Write a job file conforming to this cluster.
                          Options: ['pbs']""",
                          default='pbs')
    config_p.add_argument('--parallel', '-p',
                          help="""How to do parallel.
                          Options: 'parallel': use GNU parallel
                                 'serial': Do things serially""",
                          default='parallel')
    config_p.add_argument('--out_fn', '-o',
                          help="""Prefix for files to write.""",
                          default='maccel')
    config_p.add_argument('--n_copy', '-n',
                          help="""Number of copies to run.""",
                          type=int,
                          default=1)

    config_sp = config_p.add_subparsers(dest='config_type')
    config_sp.required = True

    config_ala = config_sp.add_parser('alanine')
    config_ala.set_defaults(config='alanine')

    args = parser.parse_args()

    print('--------------')
    for k, v in vars(args).items():
        print("{}:\t{}".format(k, v))
    print('--------------\n')

    args.func(args)


if __name__ == "__main__":
    parse()