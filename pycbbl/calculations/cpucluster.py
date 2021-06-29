def singleJob(job, script_name=None, job_name=None, cpus=1, output=None):

    if job_name == None:
        raise ValueError('job_name == None. You need to specify a name for the job')
    if output == None:
        output = job_name
    if script_name == None:
        script_name = 'slurm_submission.sh'

    #Write slurm script
    with open(script_name, 'w') as sf:
        sf.write('#!/bin/bash\n')
        sf.write('#SBATCH --job-name='+job_name+'\n')
        sf.write('#SBATCH -n '+str(cpus)+'\n')
        sf.write('#SBATCH --output='+output+'_%a_%A.out\n')
        sf.write('#SBATCH --error='+output+'_%a_%A.err\n')
        sf.write('\n')

        sf.write(job+'\n')

def jobArrays(jobs, script_name=None, job_name=None, cpus=1, gpu=False, output=None):

    if job_name == None:
        raise ValueError('job_name == None. You need to specify a name for the job')
    if output == None:
        output = job_name
    if script_name == None:
        script_name = 'slurm_array.sh'

    #Write jobs as array
    with open(script_name,'w') as sf:
        sf.write('#!/bin/bash\n')
        sf.write('#SBATCH --job-name='+job_name+'\n')
        sf.write('#SBATCH --nodes=1\n')
        sf.write('#SBATCH --ntasks=1\n')
        sf.write('#SBATCH --cpus-per-task='+str(cpus)+'\n')
        if gpu:
            sf.write('#SBATCH --gres=gpu:1'+'\n')
        sf.write('#SBATCH --array=1-'+str(len(jobs))+'\n')
        sf.write('#SBATCH --output='+output+'_%a_%A.out\n')
        sf.write('#SBATCH --error='+output+'_%a_%A.err\n')
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
