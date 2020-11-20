def jobArrays(jobs, script_name=None, job_name=None, cpus=1, gpu=False,
              threads=None, output=None, mail=None, time=72, modules=None,
              conda_env=None, unload_modules=None):

    if job_name == None:
        raise ValueError('job_name == None. You need to specify a name for the job')
    if output == None:
        output = job_name
    if script_name == None:
        script_name = 'slurm_array.sh'
    if mail == None:
        mail = 'martinfloor@gmail.com'
    if modules != None:
        if isinstance(modules, str):
            modules = [modules]
        if not isinstance(modules, list):
            raise ValueError('Modules to load must be given as a list or as a string (for loading one module only)')
    if unload_modules != None:
        if isinstance(unload_modules, str):
            unload_modules = [unload_modules]
        if not isinstance(unload_modules, list):
            raise ValueError('Modules to unload must be given as a list or as a string (for unloading one module only)')
    if conda_env != None:
        if not isinstance(conda_env, str):
            raise ValueError('The conda environment must be given as a string')

    #Write jobs as array
    with open(script_name,'w') as sf:
        sf.write('#!/bin/bash\n')
        sf.write('#SBATCH --job-name='+job_name+'\n')
        sf.write('#SBATCH --time='+str(time)+':00:00\n')
        sf.write('#SBATCH -n '+str(cpus)+'\n')
        if threads != None:
            sf.write('#SBATCH -c '+str(threads)+'\n')
        if gpu:
            sf.write('#SBATCH --gres=gpu:1'+'\n')
        sf.write('#SBATCH --array=1-'+str(len(jobs))+'\n')
        sf.write('#SBATCH --output='+output+'_%a_%A.out\n')
        sf.write('#SBATCH --error='+output+'_%a_%A.err\n')
        sf.write('#SBATCH --mail-user='+mail+'\n')
        sf.write('#SBATCH --mail-type=END,FAIL\n')
        sf.write('\n')

        if unload_modules != None:
            for module in unload_modules:
                sf.write('module unload '+module+'\n')
            sf.write('\n')
        if modules != None:
            for module in modules:
                sf.write('module load '+module+'\n')
            sf.write('\n')
        if conda_env != None:
            sf.write('conda activate '+conda_env+'\n')
            sf.write('\n')

    for i in range(len(jobs)):
        with open(script_name,'a') as sf:
            sf.write('if [[ $SLURM_ARRAY_TASK_ID = '+str(i+1)+' ]]; then\n')
            sf.write(jobs[i])
            if jobs[i].endswith('\n'):
                sf.write('fi\n')
            else:
                sf.write('\nfi\n')
            sf.write('\n')

    if conda_env != None:
        with open(script_name,'a') as sf:
            sf.write('conda deactivate \n')
            sf.write('\n')

def singleJob(job, script_name=None, job_name=None, partition='class_a', cpus=24, time=120,
              output=None, mail=None, modules=None, conda_env=None):

    available_partitions = ['class_a']

    if job_name == None:
        raise ValueError('job_name == None. You need to specify a name for the job')
    if output == None:
        output = job_name
    if script_name == None:
        script_name = 'sub_script.sh'
    if mail == None:
        mail = 'martinfloor@gmail.com'
    if modules != None:
        if isinstance(modules, str):
            modules = [modules]
        if not isinstance(modules, list):
            raise ValueError('Modules to load must be given as a list or as a string (for loading one module only)')
    if conda_env != None:
        if not isinstance(conda_env, str):
            raise ValueError('The conda environment must be given as a string')

    #Write slurm script
    with open(script_name,'w') as sf:
        sf.write('#!/bin/bash\n')
        sf.write('#SBATCH --job-name='+job_name+'\n')
        sf.write('#SBATCH --partition='+partition+'\n')
        sf.write('#SBATCH --time='+str(time)+':00:00\n')
        sf.write('#SBATCH -n '+str(cpus)+'\n')
        sf.write('#SBATCH --output='+output+'_%a_%A.out\n')
        sf.write('#SBATCH --error='+output+'_%a_%A.err\n')
        sf.write('#SBATCH --mail-user='+mail+'\n')
        sf.write('#SBATCH --mail-type=END,FAIL\n')
        sf.write('\n')

        if modules != None:
            for module in modules:
                sf.write('module load '+module+'\n')
            sf.write('\n')

        if conda_env != None:
            sf.write('conda activate '+conda_env+'\n')
            sf.write('\n')

        sf.write(job+'\n')

        if conda_env != None:
            sf.write('conda deactivate \n')
            sf.write('\n')
