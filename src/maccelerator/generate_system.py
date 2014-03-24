'''
Created on Mar 23, 2014

@author: harrigan

Code for setting up a directory structure and jobs for an accelerator run
'''

import os
import shutil
import stat
import logging as log

PBS_HEADER = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime={hours}:00:00
#PBS -l mem=8gb
#PBS -j oe
#PBS -M harrigan@stanford.edu
#PBS -m ae


cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

"""

SIMULATE_JOB = PBS_HEADER + """
for i in {{{start_i}..{end_i}}}
do
    maccel.py run --round {round_i} --traj $i --n_spt {n_spt} --report {report} &> jobs/{job_fn}.log
done
"""

MODEL_JOB = PBS_HEADER + """
maccel.py model --round {round_i} --lagtime {lagtime} --n_tpr {n_tpr} &> jobs/{job_fn}.log
"""

SIMULATE_SUBMIT = """S{traj_i}=`mqsub {dep} jobs/{job_fn}.job`"""

MODEL_SUBMIT = """M{round_i}=`mqsub {dep} jobs/{job_fn}.job`"""


def create_file_structure(proj_dir, seed_structures):
    """Create folders and files to initialize a directory for a run.

    :proj_dir: project directory name
    :seed_structures: path to initial seed structures, which will be copied in

    """

    # Make project
    log.info('Creating directory %s', proj_dir)
    os.mkdir(proj_dir)

    # Make starting structures
    os.mkdir(os.path.join(proj_dir, 'sstates'))
    shutil.copy(
        seed_structures,
        os.path.join(
            proj_dir,
            'sstates',
            'round-0.h5'))

    # Make jobs and trajs dir
    os.mkdir(os.path.join(proj_dir, 'jobs'))
    os.mkdir(os.path.join(proj_dir, 'trajs'))


def write_dep_jobs(proj_dir, args):
    """Write files for a dependency run.

    :proj_dir: project directory
    :start_from: round from which to start, none for starting from beginning
    :args: from command line
    """

    start_from = args.start_from
    if start_from is None:
        create_file_structure(proj_dir, args.seed_structures)
        start_from = 0
    else:
        # We used to try to recover the qsub id so that we could restart
        # before it was finished. But if the job no longer exists, the jobs
        # never start
        pass

    submit_lines = []
    # Iterate over rounds
    for round_i in range(start_from, args.n_round):

        # Make the trajectory subdirectory
        try:
            os.mkdir(os.path.join(proj_dir, 'trajs', 'round-%d' % round_i))
        except OSError:
            pass

        # Initialize dependency spec for modeling
        m_dep = '-W depend=afterok'

        # Make simulate jobs
        for traj_i in range(0, args.n_tpr, args.n_tpj):
            job_fn = 'round-{round_i}_traj-{traj_i}'.format(round_i=round_i,
                                                            traj_i=traj_i)

            with open(os.path.join(proj_dir, 'jobs', "%s.job" % job_fn), 'w') as job_f:
                job_f.write(SIMULATE_JOB.format(round_i=round_i,
                                                n_spt=args.n_spt * args.report,
                                                report=args.report, job_fn=job_fn,
                                                hours=1,
                                                start_i=traj_i,
                                                end_i=traj_i + args.n_tpj - 1))

                if round_i > start_from:
                    dep = '-W depend=afterok:$M{pti}'.format(pti=round_i - 1)
                else:
                    dep = ''

                submit_lines += [SIMULATE_SUBMIT.format(traj_i=traj_i,
                                                        job_fn=job_fn,
                                                        dep=dep)]
                m_dep += ':$S{traj_i}'.format(traj_i=traj_i)

        # Make model job
        job_fn = 'round-{round_i}_model'.format(round_i=round_i)
        with open(os.path.join(proj_dir, 'jobs', "%s.job" % job_fn), 'w') as job_f:
            job_f.write(MODEL_JOB.format(round_i=round_i,
                                         lagtime=args.lagtime,
                                         n_tpr=args.n_tpr, hours=3,
                                         job_fn=job_fn))
            submit_lines += [MODEL_SUBMIT.format(job_fn=job_fn,
                                                 dep=m_dep,
                                                 round_i=round_i)]

    with open(os.path.join(proj_dir, 'submit.sh'), 'w') as sub_f:
        sub_f.write('\n'.join(submit_lines))

    # Make executable
    st = os.stat(os.path.join(proj_dir, 'submit.sh'))
    os.chmod(os.path.join(proj_dir, 'submit.sh'), st.st_mode | stat.S_IEXEC)


def write_one_job(proj_dir, args):
    """Write files where there will only be one job."""
    # TODO
    
    create_file_structure(proj_dir, args.seed_structures)
    




