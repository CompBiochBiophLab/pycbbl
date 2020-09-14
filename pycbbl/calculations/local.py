import os

def parallel(jobs, cpus=6, script_name='commands'):
    """
    Generates scripts to run jobs simultaneously in N Cpus in a local computer,
    i.e. without a job manager. The input jobs must be a list representing each
    job to execute as a string.

    Two different scripts are written to execute the jobs in bash language. For
    example, if the script_name variable is set to commands and the cpus to 4, five
    scripts will be written:

    - commands
    - commands_0
    - commands_1
    - commands_2
    - commands_3

    The jobs to execute are distributed into the numbered scripts. Each numbered
    script contains a sub set of jobs that will be executed in a sequential manner.
    The numberless script execute all the numbered scripts in the background, using
    the nohup command, and redirecting the output to different files for each numbered
    script. To execute the jobs is only necessary to execute:

    'bash commands'

    Parameters
    ----------
    jobs : list
        List of strings containing the commands to execute jobs.
    cpus : int
        Number of CPUs to use in the execution.
    script_name : str
        Name of the output scripts to execute the jobs.
    """
    # Write parallel execution scheme #

    dJobs = int(len(jobs)/(cpus))
    rJobs = int(len(jobs)%(cpus))

    zf = len(str(cpus))
    count = 0
    for i in range(cpus):
        with open(script_name+'_'+str(i).zfill(zf),'w') as sf:
            for j in range(dJobs):
                sf.write(jobs[count])
                count += 1

    for i in range(cpus):
        with open(script_name+'_'+str(i).zfill(zf),'a') as sf:
            for j in range(rJobs):
                sf.write(jobs[count])
                count += 1

    with open(script_name,'w') as sf:
        sf.write('for script in '+script_name+'_'+'?'*zf+'; do nohup bash $script &> ${script%.*}.nohup& done\n')

def multipleGPUSimulations(jobs, parallel=3, gpus=4, script_name='gpu_commands'):
    """
    Generates scripts to run jobs simultaneously in N GPUs and X cpus (parallel option)
    in a local computer, i.e. without a job manager. The input jobs must be a list
    representing each job to execute as a string. Each job string must contain the
    substring 'GPUID' which will be replaced by an integer representing the GPU
    to use. When more than one simultaneous jobs are to be executed in the same
    GPU, the 'parallel' option will distribute the jobs into N cpus more giving
    a total of 'N*parallel' executions, with a 'parallel' number of jobs being in
    executed in a single GPU.

    Two different scripts are written to execute the jobs in bash language. For
    example, if the script_name variable is set to commands and the cpus to 4, five
    scripts will be written:

    - commands
    - commands_0
    - commands_1
    - commands_2
    - commands_3

    The jobs to execute are distributed into the numbered scripts. Each numbered
    script contains a sub set of jobs that will be executed in a sequential manner.
    The numberless script execute all the numbered scripts in the background, using
    the nohup command, and redirecting the output to different files for each numbered
    script. To execute the jobs is only necessary to execute:

    'bash commands'

    Parameters
    ----------
    jobs : list
        List of strings containing the commands to execute jobs.
    gpus : int
        Number of GPUs to use in the execution.
    parallel : int
        Number of parallel executions into the same GPU.
    script_name : str
        Name of the output scripts to execute the jobs.
    """
    # Write parallel execution scheme #
    dJobs = int(len(jobs)/(gpus*parallel))
    rJobs = int(len(jobs)%(gpus*parallel))
    gpus_count = 0
    count = 0

    if len(jobs) <= gpus*parallel:
        zf = len(str(len(jobs)))
    else:
        zf = len(str(gpus*parallel))

    for i in range(gpus*parallel):
        gpuId = gpus_count
        with open(script_name+'_'+str(i).zfill(zf), 'w') as sf:
            for j in range(dJobs):
                sf.write(jobs[count].replace('GPUID',str(gpuId)))
                count += 1
        gpus_count += 1
        if gpus_count == gpus:
            gpus_count = 0


    for i in range(rJobs):
        gpuId = gpus_count
        with open(script_name+'_'+str(i).zfill(zf),'a') as sf:
            sf.write(jobs[count].replace('GPUID',str(gpuId)))
        count += 1
        gpus_count += 1
        if gpus_count == gpus:
            gpus_count = 0

    if dJobs == 0:
        for i in range(rJobs, gpus*parallel):
            os.remove(script_name+'_'+str(i).zfill(zf))

    with open(script_name,'w') as sf:
        sf.write('for script in '+script_name+'_'+'?'*zf+'; do nohup bash $script &> ${script%.*}.nohup& done\n')